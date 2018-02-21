\name{plot.flexgam}
\alias{plot.flexgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots object of class flexgam.
}
\description{
Plots resulting response (link) function or estimated smooth effects of objects of class \code{flexgam}.
}
\usage{
\method{plot}{flexgam}(x, type=c("response","covariate"), ci = TRUE,  rug = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class flexgam.
}
  \item{type}{
Whether the response function or the smooth covariate effects should be plotted.
}
  \item{ci}{
Include confidence intervals?
}
  \item{rug}{
Include rug plots?
}
  \item{...}{
Standard plot add-ons should work. ylab and xlab check for the length of the smooth predictor.
}
}
\details{
Plots either the estimated response function or the estimated smooth covariate effects. Valid confidence intervals are currently only plotted for P-splines.
}
\value{
Plots
}
\references{
Spiegel, Elmar, Thomas Kneib and Fabian Otto-Sobotka. Generalized additive models with flexible response functions. Statistics and Computing (2017). https://doi.org/10.1007/s11222-017-9799-6
}
\author{
Elmar Spiegel
}
\note{
If other smooth effects than P-splines are used in the formula the plot function internally calls \code{mgcv::plot.gam} to show their effect. Then P-splines are plotted twice.
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

plot(model_2,type="response")
plot(model_2,type="covariate")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
