\name{response.flexgam}
\alias{response.flexgam}
\alias{response}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicts values for the object of class flexgam
}
\description{
Predicts values for the given object, if a new predictor is given.
}
\usage{
\method{response}{flexgam}(object, linear.predictor = NULL, ...)
response(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class flexgam.
}
  \item{linear.predictor}{
New predictor, which has to be centered and scaled properly.
}
  \item{...}{
Currently not used
}
}
\details{
Calculates predicted values for the given model. Used to get the estimated response function separately from the covariate effects. For the same predictor \code{predict.flexgam(\ldots,type="response")} and \code{response(\ldots)} are equal.
}
\value{
Numeric vector of fitted values.
}
\references{
Spiegel, Elmar, Thomas Kneib and Fabian Otto-Sobotka. Generalized additive models with flexible response functions. Statistics and Computing (2017). https://doi.org/10.1007/s11222-017-9799-6
}
\author{
Elmar Spiegel
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{flexgam}}, \code{\link{deviance.flexgam}}, \code{\link{predict.flexgam}}
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


set.seed(2)
n <- 1000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
eta_orig <- -1 + 2*sin(6*x1) + exp(x2) + x3
pi_orig <- pgamma(eta_orig, shape=2, rate=sqrt(2))
y <- rbinom(n,size=1,prob=pi_orig)

newData <- data.frame(y,x1,x2,x3)



fitted_1 <- predict(model_2, newdata=newData, type="response")
predictor1 <- predict(model_2, newdata=newData, type="linear.predictor")
fitted_2 <- response(model_2, linear.predictor=predictor1)
all.equal(fitted_1,fitted_2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
