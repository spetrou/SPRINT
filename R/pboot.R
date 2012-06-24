##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

sample0 <- function(x, ...) x[sample.int(length(x), ...)]
bsample <- function(x, ...) x[sample.int(length(x), replace = TRUE, ...)]

isMatrix <- function(x) length(dim(x)) == 2L

pboot <- function(data, statistic, R, sim = "ordinary",
                  stype = c("i", "f", "w"),
                  strata  =  rep(1, n), L = NULL, m = 0, weights = NULL,
                  ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ...)
{  
#
# R replicates of bootstrap applied to  statistic(data)
# Possible sim values are: "ordinary", "balanced", "antithetic",
#                     "permutation", "parametric"
# Various auxilliary functions find the indices to be used for the
# bootstrap replicates and then this function loops over those replicates.
#
  call <- match.call()
  stype <- match.arg(stype)
  
  if (simple && (sim != "ordinary" || stype != "i" || sum(m))) {
    warning("'simple=TRUE' is only valid for 'sim=\"ordinary\", stype=\"i\", n=0, so ignored")
    simple <- FALSE
  }
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    n <- NROW(data)
  if ((n == 0) || is.null(n))
    stop("no data in call to boot")
  temp.str <- strata
  strata <- tapply(seq_len(n),as.numeric(strata))
  t0 <- if (sim != "parametric") {
    if ((sim == "antithetic") && is.null(L))
      L <- empinf(data = data, statistic = statistic,
                  stype = stype, strata = strata, ...)
    if (sim != "ordinary") m <- 0
    else if (any(m < 0)) stop("negative value of m supplied")
    if ((length(m) != 1L) && (length(m) != length(table(strata))))
      stop("length of m incompatible with strata")
    if ((sim == "ordinary") || (sim == "balanced")) {
      if (isMatrix(weights) && (nrow(weights) != length(R)))
        stop("dimensions of R and weights do not match")}
    else weights <- NULL
    if (!is.null(weights))
      weights <- t(apply(matrix(weights, n, length(R), byrow = TRUE),
                         2L, normalize, strata))
    if (!simple) i <- index.array(n, R, sim, strata, m, L, weights)
    
    original <- if (stype == "f") rep(1, n)
    else if (stype == "w") {
      ns <- tabulate(strata)[strata]
      1/ns
    } else seq_len(n)
    
    t0 <- if (sum(m) > 0L) statistic(data, original, rep(1, sum(m)), ...)
    else statistic(data, original, ...)
    rm(original)
        t0
  } else # "parametric"
  statistic(data, ...)
  
  pred.i <- NULL
  fn <- if (sim == "parametric") {
    ## force promises, so values get sent by parallel
    ran.gen; data; mle
    function(r) {
      dd <- ran.gen(data, mle)
      statistic(dd, ...)
    }
  } else {
    if (!simple && ncol(i) > n) {
      pred.i <- as.matrix(i[ , (n+1L):ncol(i)])
      i <- i[, seq_len(n)]
    }
    if (stype %in% c("f", "w")) {
      f <- freq.array(i)
      rm(i)
      if (stype == "w") f <- f/ns
      if (sum(m) == 0L) function(r) statistic(data, f[r,  ], ...)
      else function(r) statistic(data, f[r, ], pred.i[r, ], ...)
    } else if (sum(m) > 0L)
      function(r) statistic(data, i[r, ], pred.i[r,], ...)
    else if (simple)
      function(r)
        statistic(data,
                  index.array(n, 1, sim, strata, m, L, weights), ...)
    else function(r) statistic(data, i[r, ], ...)
  }
  RR <- sum(R)
  
  res <- .Call("pboot",
               as.list(seq_len(RR)),
               fn)

  t.star <- matrix(, RR, length(t0))
  for(r in seq_len(RR)) t.star[r, ] <- res[[r]]
  
  if (is.null(weights)) weights <- 1/tabulate(strata)[strata]
  boot.return(sim, t0, t.star, temp.str, R, data, statistic, stype, call,
              seed, L, m, pred.i, weights, ran.gen, mle)
}

boot.return <- function(sim, t0, t, strata, R, data, stat, stype, call,
			seed, L, m, pred.i, weights, ran.gen, mle)
                                        #
# Return the results of a bootstrap in the form of an object of class
# "boot".
#
{
  out <- list(t0=t0, t=t, R=R, data=data, seed=seed,
              statistic=stat, sim=sim, call=call)
  if (sim == "parametric")
    out <- c(out, list(ran.gen=ran.gen, mle=mle))
  else if (sim == "antithetic")
    out <- c(out, list(stype=stype, strata=strata, L=L))
  else if (sim == "ordinary") {
    if (sum(m) > 0)
      out <- c(out, list(stype=stype, strata=strata,
                         weights=weights, pred.i=pred.i))
    else 	out <- c(out, list(stype=stype, strata=strata,
                                   weights=weights))
  } else if (sim == "balanced")
    out <- c(out, list(stype=stype, strata=strata,
                           weights=weights ))
  else
    out <- c(out, list(stype=stype, strata=strata))
  class(out) <- "boot"
  out
}

index.array <- function(n, R, sim, strata=rep(1,n), m=0, L=NULL, weights=NULL)
{
#
#  Driver function for generating a bootstrap index array.  This function
#  simply determines the type of sampling required and calls the appropriate
#  function.
#
  indices <- NULL
  if (is.null (weights)) {
    if (sim == "ordinary") {
      indices <- ordinary.array(n, R, strata)
      if (sum(m) > 0)
        indices <- cbind(indices, extra.array(n, R, m, strata))
    }
    else if (sim == "balanced")
      indices <- balanced.array(n, R, strata)
    else if (sim == "antithetic")
      indices <- antithetic.array(n, R, L, strata)
    else if (sim == "permutation")
      indices <- permutation.array(n, R, strata)
  } else {
    if (sim == "ordinary")
      indices <- importance.array(n, R, weights, strata)
    else if (sim == "balanced")
      indices <- importance.array.bal(n, R, weights, strata)
  }
  indices
}

ordinary.array <- function(n, R, strata)
{
#
# R x n array of bootstrap indices, resampled within strata.
# This is the function which generates a regular bootstrap array
# using equal weights within each stratum.
#
  inds <- as.integer(names(table(strata)))
  if (length(inds) == 1L) {
    output <- sample.int(n, n*R, replace=TRUE)
    dim(output) <- c(R, n)
  } else {
    output <- matrix(as.integer(0L), R, n)
    for(is in inds) {
      gp <- seq_len(n)[strata == is]
      output[, gp] <- if (length(gp) == 1) rep(gp, R) else bsample(gp, R*length(gp))
    }
  }
  output
}


freq.array <- function(i.array)
{
#
# converts R x n array of bootstrap indices into
# R X n array of bootstrap frequencies
#
  result <- NULL
  n <- ncol(i.array)
  result <- t(apply(i.array, 1, tabulate, n))
  result
}

permutation.array <- function(n, R, strata)
{
#
# R x n array of bootstrap indices, permuted within strata.
# This is similar to ordinary array except that resampling is
# done without replacement in each row.
#
  output <- matrix(rep(seq_len(n), R), n, R)
  inds <- as.integer(names(table(strata)))
  for(is in inds) {
    group <- seq_len(n)[strata == is]
    if (length(group) > 1L) {
      g <- apply(output[group,  ], 2L, rperm)
      output[group,  ] <- g
    }
  }
  t(output)
}

normalize <- function(wts, strata)
{
#
# Normalize a vector of weights to sum to 1 within each strata.
#
  n <- length(strata)
  out <- wts
  inds <- as.integer(names(table(strata)))
  for (is in inds) {
    gp <- seq_len(n)[strata == is]
    out[gp] <- wts[gp]/sum(wts[gp]) }
  out
}
