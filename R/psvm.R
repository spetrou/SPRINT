psvm <- function(x,
    y           = NULL,
    type        = 'C',
    kernel      = "radial",
    degree      = 3,
    gamma       = if (is.vector(x)) 1 else 1 / ncol(x),
    coef0       = 0,
    cost        = 1,
    nu          = 0.5,
    class.weights = NULL,
    cachesize   = 40.0,
    tolerance   = 0.001,
    epsilon     = 0.1,
    shrinking   = TRUE,
    cross       = 0,
    probability = FALSE,
    fitted      = TRUE,
    seed        = 1L,
    scale       = TRUE,
    na.action   = na.omit)

{
    if(inherits(x, "Matrix")) {
       library("SparseM")
	library("Matrix")
	x <- as(x, "matrix.csr")
    }
    if(inherits(x, "simple_triplet_matrix")) {
	library("SparseM")
	ind <- order(x$i, x$j)
	x <- new("matrix.csr",
	    ra = x$v[ind],
	    ja = x$j[ind],
	    ia = as.integer(cumsum(c(1, tabulate(x$i[ind])))),
	    dimension = c(x$nrow, x$ncol))
     }
     if (sparse <- inherits(x, "matrix.csr")){
	 library("SparseM")
     }
   ## determine model type
    
    type <- pmatch(type, c("C-classification",
	    "nu-classification",
	    "one-classification",
	    "eps-regression",
	    "nu-regression"), 99) - 1
    
#    cat("Type = \"", type, "\"\n" )
    
    kernel <- pmatch(kernel, c("linear",
	    "polynomial",
	    "radial",
	    "sigmoid"), 99) - 1
    
 #   cat("Kernel = \"", kernel,"\"\n" )
    
  
    formula <- inherits(x, "svm.formula")
    nac <- attr(x, "na.action")
    
 ## further parameter checks
    nr <- nrow(x)
    ly <- length(y)

    if (cross > nr)
        stop(sQuote("cross"), " cannot exceed the number of observations!")
    if (!is.vector(y) && !is.factor (y) && type != 2)
        stop("y must be a vector or a factor.")
    if (type != 2 && ly != nr)
        stop("x and y don't match.")
    if (type > 2 && !is.numeric(y))
        stop("Need numeric dependent variable for regression.")

    ## scaling, subsetting, and NA handling
    if (sparse) {

	scale <- rep(FALSE, ncol(x))
	if(!is.null(y)) na.fail(y)
	x <- SparseM::t(SparseM::t(x)) ## make shure that col-indices are sorted
		
     } else {
	x <- as.matrix(x)

	## subsetting and na-handling for matrices
        if (!formula) {
	    if (is.null(y))
                x <- na.action(x)
            else {

		#generates a frame of the input data
                df <- na.action(data.frame(y, x))
                y <- df[,1]
		x <- as.matrix(df[,-1])
	        nac <-
                    attr(x, "na.action") <-
                        attr(y, "na.action") <-
                            attr(df, "na.action")		
            }
        }

     }   
weightlabels <- 0
       	
## in case of classification: transform factors into integers
    if (type == 2) # one class classification --> set dummy
        y <- rep(1, nr)
    else
        if (is.factor(y)) {
            lev <- levels(y)
            y <- as.integer(y)
            if (!is.null(class.weights)) {
                if (is.null(names(class.weights)))
                    stop("Weights have to be specified along with their according level names !")
                weightlabels <- match (names(class.weights), lev)
                if (any(is.na(weightlabels)))
                    stop("At least one level name is missing or misspelled.")
            }
        } else {
            if (type < 3) {
                if(any(as.integer(y) != y))
                    stop("dependent variable has to be of factor or integer type for classification mode.")

		y <- as.factor(y)
                lev <- levels(y)
                y <- as.integer(y)

            } else 
		lev <- unique(y)
        }
   
    
    #nclass <- 2
    if (type < 2) nclass <- length(lev)

    if (type > 1 && length(class.weights) > 0) {
        class.weights <- NULL
        warning(sQuote("class.weights"), " are set to NULL for regression mode. For classification, use a _factor_ for ", sQuote("y"),
", or specify the correct ", sQuote("type"), " argument.")
    }

    err <- empty_string <- paste(rep(" ", 255), collapse = "")
    
#   spRows <- if (sparse) x@ia else 0
#   spCols <- if (sparse) x@ja else 0
    
    spRows <- 10
    spCols <- 10
    ret <- .Call("psvm", 
        
        ## data
	as.double  (t(x)),
	as.integer (nr), 
	as.integer(ncol(x)),
	as.double  (y),
	as.integer (nclass),
	as.integer (cross),
	
        ## sparse index info
	as.integer (spRows),
	as.integer (spCols),

	## parameters
	as.integer (type),
	as.integer (kernel),
	as.integer (degree),
	as.double  (gamma),
	as.double  (coef0),
	as.double  (cost),
	as.double  (nu),
	as.integer (weightlabels),
	as.double  (class.weights),
	as.integer (length (class.weights)),
	as.double  (cachesize),
	as.double  (tolerance),
	as.double  (epsilon),
	as.integer (shrinking),
	as.integer (sparse),
	as.integer (probability),
	as.integer (seed)

        ## results

#	nclasses = integer  (1),
#	nr       = integer  (1), # nr of support vectors
#	index    = integer  (nr),
#	labels   = integer  (nclass),
#	nSV      = integer  (nclass),

#	rho      = double   (nclass * (nclass - 1) / 2),
#	coefs    = double   (nr * (nclass - 1)),
#	sigma    = double   (1),
#	probA    = double   (nclass * (nclass - 1) / 2),
#	probB    = double   (nclass * (nclass - 1) / 2),
#	
#	cresults = double   (cross),
#	ctotal1  = double   (1),
#	ctotal2  = double   (1),
#	error    = err
    )
   
    return(ret)
}
