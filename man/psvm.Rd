\name{psvm}
\alias{psvm}
\title{
Parallel Support Vector Machines
}
\description{
 A parallel implementation of the svm function from R package e1071. 
 \code{psvm} is used to train a support vector machine. It can be used to carry
out general regression and classification (of nu and epsilon-type), as
well as density-estimation. A formula interface is provided.
}
\usage{
psvm(x, y = NULL, type = "C", kernel = "radial", degree = 3, 
gamma = if (is.vector(x)) 1 else 1/ncol(x), coef0 = 0, cost = 1, 
nu = 0.5, class.weights = NULL, cachesize = 40, tolerance = 0.001, 
epsilon = 0.1, shrinking = TRUE, cross = 0, probability = FALSE, 
fitted = TRUE, seed = 1L, scale = TRUE, na.action = na.omit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%  \item{formula}{a symbolic description of the model to be fit.}
%  \item{data}{an optional data frame containing the variables in the model.
%          By default the variables are taken from the environment which
%          \sQuote{psvm} is called from.}
  \item{x}{a data matrix, a vector, or a sparse matrix (object of class
    \code{\link[Matrix]{Matrix}} provided by the \pkg{Matrix} package,
    or of class \code{\link[SparseM]{matrix.csr}}
    provided by the \pkg{SparseM} package, or of class
    \code{\link[slam]{simple_triplet_matrix}} provided by the \pkg{slam}
    package).}
  \item{y}{a response vector with one label for each row/component of
    \code{x}. Can be either a factor (for classification tasks)
    or a numeric vector (for regression).}
  \item{scale}{A logical vector indicating the variables to be
    scaled. If \code{scale} is of length 1, the value is recycled as
    many times as needed.
    Per default, data are scaled internally (both \code{x} and \code{y}
    variables) to zero mean and unit variance. The center and scale
    values are returned and used for later predictions.}
  \item{type}{\code{psvm} can be used as a classification
    machine, as a regression machine, or for novelty detection.
    Depending of whether \code{y} is
    a factor or not, the default setting for \code{type} is \code{C-classification} or \code{eps-regression}, respectively, but may be overwritten by setting an explicit value.\cr
    Valid options are:
    \itemize{
      \item \code{C-classification}
      \item \code{nu-classification}
      \item \code{one-classification} (for novelty detection)
      \item \code{eps-regression}
      \item \code{nu-regression}
    }
  }
  \item{kernel}{the kernel used in training and predicting. You
    might consider changing some of the following parameters, depending
    on the kernel type.\cr
    \describe{
      \item{linear:}{\eqn{u'v}{u'*v}}
      \item{polynomial:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
      \item{radial basis:}{\eqn{e^(-\gamma |u-v|^2)}{exp(-gamma*|u-v|^2)}}
      \item{sigmoid:}{\eqn{tanh(\gamma u'v + coef0)}{tanh(gamma*u'*v + coef0)}}
      }
    }
  \item{degree}{parameter needed for kernel of type \code{polynomial} (default: 3)}
  \item{gamma}{parameter needed for all kernels except \code{linear}
    (default: 1/(data dimension))}
  \item{coef0}{parameter needed for kernels of type \code{polynomial}
    and \code{sigmoid} (default: 0)}
  \item{cost}{cost of constraints violation (default: 1)---it is the
    \sQuote{C}-constant of the regularization term in the Lagrange formulation.}
  \item{nu}{parameter needed for \code{nu-classification},
    \code{nu-regression}, and \code{one-classification}}
  \item{class.weights}{a named vector of weights for the different
    classes, used for asymmetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.}
  \item{cachesize}{cache memory in MB (default 40)}
  \item{tolerance}{tolerance of termination criterion (default: 0.001)}
  \item{epsilon}{epsilon in the insensitive-loss function (default: 0.1)}
  \item{shrinking}{option whether to use the shrinking-heuristics
    (default: \code{TRUE})}
  \item{cross}{if a integer value k>0 is specified, a k-fold cross
    validation on the training data is performed to assess the quality
    of the model: the accuracy rate for classification and the Mean
    Squared Error for regression}
  \item{fitted}{logical indicating whether the fitted values should be computed
    and included in the model or not (default: \code{TRUE})}
  \item{probability}{logical indicating whether the model should
    allow for probability predictions.}
  \item{seed}{integer seed for libsvm (used for cross-validation and
    probability prediction models).}
  \item{\dots}{additional parameters for the low level fitting function
    \code{psvm.default}}
%  \item{subset}{An index vector specifying the cases to be used in the
%          training sample.  (NOTE: If given, this argument must be
%          named.)}
  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.omit}, which leads to rejection of cases
          with missing values on any required variable. An alternative
	  is \code{na.fail}, which causes an error if \code{NA} cases
	  are found. (NOTE: If given, this argument must be named.)}	
}
\value{
  An object of class \code{"psvm"} containing the fitted model, including:
  \item{SV}{The resulting support vectors (possibly scaled).}
  \item{index}{The index of the resulting support vectors in the data
    matrix. Note that this index refers to the preprocessed data (after
    the possible effect of \code{na.omit} and \code{subset})}
  \item{coefs}{The corresponding coefficients times the training labels.}
  \item{rho}{The negative intercept.}
  \item{sigma}{In case of a probabilistic regression model, the scale
    parameter of the hypothesized (zero-mean) laplace distribution estimated by
    maximum likelihood.}
  \item{probA, probB}{numeric vectors of length k(k-1)/2, k number of
    classes, containing the parameters of the logistic distributions fitted to
    the decision values of the binary classifiers (1 / (1 + exp(a x + b))).}
}
\details{
  For multiclass-classification with k levels, k>2, \code{libsvm} uses the
  \sQuote{one-against-one}-approach, in which k(k-1)/2 binary classifiers are
  trained; the appropriate class is found by a voting scheme.
  
  \code{libsvm} internally uses a sparse data representation, which is 
  also high-level supported by the package \pkg{SparseM}.
  
  If the predictor variables include factors, the formula interface must be used to get a
  correct model matrix.

  \code{plot.svm} allows a simple graphical
  visualization of classification models.

  The probability model for classification fits a logistic distribution
  using maximum likelihood to the decision values of all binary
  classifiers, and computes the a-posteriori class probabilities for the
  multi-class problem using quadratic optimization. The probabilistic
  regression model assumes (zero-mean) laplace-distributed errors for the
  predictions, and estimates the scale parameter using maximum likelihood.
}


\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    Based on svm by David Meyer (based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin)
    }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{SPRINT}}
  \code{\link{svm}}
  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}
\keyword{utilities}
\keyword{interface}