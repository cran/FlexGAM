\name{flexgam_outputs}
\alias{print.flexgam}
\alias{summary.flexgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints object of class flexgam
}
\description{
Prints information for objects of class \code{flexgam}.
}
\usage{
\method{print}{flexgam}(x, ...)
\method{summary}{flexgam}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class \code{flexgam}.
}
  \item{x}{
Object of class \code{flexgam}.
}
  \item{...}{
Currently not used
}
}
\details{
\code{print} prints a short list of information about the given object.
\code{summary} prints the same list as \code{print} and adds the values of the estimated coefficients and the corresponding standard deviations.
}
\value{
None, only printing of object details.
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

print(model_2)
summary(model_2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
