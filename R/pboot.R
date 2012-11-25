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
      L <- boot::empinf(data = data, statistic = statistic,
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
                         2L, boot:::normalize, strata))
    if (!simple) i <- boot:::index.array(n, R, sim, strata, m, L, weights)
    
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
      f <- boot:::freq.array(i)
      rm(i)
      if (stype == "w") f <- f/ns
      if (sum(m) == 0L) function(r) statistic(data, f[r,  ], ...)
      else function(r) statistic(data, f[r, ], pred.i[r, ], ...)
    } else if (sum(m) > 0L)
      function(r) statistic(data, i[r, ], pred.i[r,], ...)
    else if (simple)
      function(r)
        statistic(data,
                  boot:::index.array(n, 1, sim, strata, m, L, weights), ...)
    else function(r) statistic(data, i[r, ], ...)
  }
  RR <- sum(R)
  
  res <- .Call("pboot",
               as.list(seq_len(RR)),
               fn)

  t.star <- matrix(, RR, length(t0))
  for(r in seq_len(RR)) t.star[r, ] <- res[[r]]

  if (is.null(weights)) weights <- 1/tabulate(strata)[strata]
  boot:::boot.return(sim, t0, t.star, temp.str, R, data, statistic, stype,
                     call, seed, L, m, pred.i, weights, ran.gen, mle)
}
