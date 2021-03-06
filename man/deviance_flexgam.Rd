\name{deviance.flexgam}
\alias{deviance.flexgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deviance of FlexGAM model
}
\description{
Calculating the (predictive) deviance of the model.
}
\usage{
\method{deviance}{flexgam}(object, newdata=NULL, use_penalty=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class \code{flexgam}.
}
  \item{newdata}{
Data to estimate the (predictive) deviance.
}
  \item{use_penalty}{
If the deviance should be penalized according to the smoothing penalties.
}
  \item{...}{
Currently not used
}
}
\details{
Calculates the model deviance of the object for the given data. To get the same value as for \code{object$deviance} you need to set \code{use_penalty = TRUE}. This is due to the fact that the deviance element of the object is the penalized deviance used for step halving.
}
\value{
Estimated deviance
}
\references{
Spiegel, Elmar, Thomas Kneib and Fabian Otto-Sobotka. Generalized additive models with flexible response functions. Statistics and Computing (2017). https://doi.org/10.1007/s11222-017-9799-6
}
\author{
Elmar Spiegel
}
\note{
To get numeric stability the fitted values are truncated (\code{min_mu_k}) to achieve numeric stability.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{flexgam}}, \code{\link{predict.flexgam}}
}
\examples{
set.seed(1)
n <- 1000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
eta_orig <- -1 + 2*sin(6*x1) + exp(x2) + x3
pi_orig <- pgamma(eta_orig, shape=2, rate=sqrt(2))
y <- rbinom(n,size=1,prob=pi_orig)

Data <- data.frame(y,x1,x2,x3)
formula <- y ~ s(x1,k=20,bs="ps") + s(x2,k=20,bs="ps") + x3

# Fix smoothing parameters to save computational time.
control2 <- list("fix_smooth" = TRUE, "quietly" = TRUE, "sm_par_vec" = 
                     c("lambda" = 100, "s(x1)" = 2000, "s(x2)" = 9000))
set.seed(2)
model_2 <- flexgam(formula=formula, data=Data, type="FlexGAM2", 
                   family=binomial(link=logit), control = control2)
(deviance           <- deviance(model_2, newdata=Data))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
