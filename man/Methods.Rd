\name{Methods}
\alias{coefficients.flexgam}
\alias{coef.flexgam}
\alias{fitted.flexgam}
\alias{fitted.values.flexgam}
\alias{residuals.flexgam}
\alias{resid.flexgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract model parameters
}
\description{
These are the standard functions to extract parameters of the estimated object.
}
\usage{
\method{coefficients}{flexgam}(object, ...)
\method{coef}{flexgam}(object, ...)
\method{fitted}{flexgam}(object, ...)
\method{fitted.values}{flexgam}(object, ...)
\method{residuals}{flexgam}(object, ...)
\method{resid}{flexgam}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class \code{flexgam}.
}
  \item{...}{
Currently not used
}
}
\details{
These functions extract the coefficients, fitted values or residuals of the given object.
}
\value{
Coefficients, fitted values or residuals of the given object.
}
\references{
Spiegel, Elmar, Thomas Kneib and Fabian Otto-Sobotka. Generalized additive models with flexible response functions. Statistics and Computing (2017). https://doi.org/10.1007/s11222-017-9799-6
}
\author{
Elmar Spiegel
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{flexgam}} 
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

coefficients(model_2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
