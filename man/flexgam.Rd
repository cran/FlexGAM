\name{flexgam}
\alias{flexgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of generalized additive model with flexible response function
}
\description{
Main function to estimate generalized additive models with flexible response function. Currently the response must be binomial, Poisson, Gaussian or Gamma distributed.
}
\usage{
flexgam(formula, data, type = c("FlexGAM1","FlexGAM2","FlexGAM1n",
        "FlexGAM2n"), family, control = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula of the covariate effects. The formula must be in the design of the \code{mgcv} package.
}
  \item{data}{
Data to fit the model.
}
  \item{type}{
Should the response function be estimated completely flexible (\code{"FlexGAM2"}) or with a surrounding canonical link (\code{"FlexGAM1"})? \code{"FlexGAM1n"} and \code{"FlexGAM2n"} are similar to \code{"FlexGAM1"} and \code{"FlexGAM2"}, but without monotonicity constraint. Therefore we do not recommend to use \code{"FlexGAM1n"} and \code{"FlexGAM2n"}.
}
  \item{family}{
Family of the data. Currently only \code{binomial(link="logit")},  \code{poisson(link="log")}, \code{gaussian(link="identity")} and \code{Gamma(link="log")} are possible.
}
  \item{control}{
Control parameters to fit the model. The default values are described in \cr \code{\link{match_flexgam_control}}.
}
}
\details{
To reduce the bias of miss-specified response functions the function estimates the response function jointly with the covariate effects. The covariate effects are build similar to the standard mgcv::gam, while the response function is either estimated as a strictly  monotone P-spline or a combination of the canonical  link and a transformation of the "linear"-predictor. In the outer loop the response function is estimated, while in the inner loop a modified version of the Fisher-Scoring algorithm is applied to get the covariate effects. In the algorithm step-halving is applied. Identifiability is achieved due to at least two smooth effects and scaling of the predictors.
}
\value{
Object of class \code{flexgam}. The list includes the 
\itemize{
\item \code{f_k}: The estimated response function.
\item \code{gam_temp}: The final step of the Fisher-Scoring algorithm, so the weighted linear model based on the mgcv-package.
\item \code{sm_par_vec}: The estimated smoothing parameters.
\item \code{coefficients}: The coefficients of the predictor as well as of the response function.
\item \code{se}: The standard deviation of the coefficients.
\item \code{mean_eta_k}, \code{sd_eta_k}: Some information about the scaling of the predictor.
%\item Information about the convergence (\code{conv_in}, \code{conv_out}, \code{conv_out_dev}, \code{iter_out}, \code{iter_in})
\item \code{control}: The applied control parameters.
\item \code{control_input}: The control parameters in the input
\item \code{details}: Information about the occurrence  of modifications due to extreme values and the convergence. As well as information of the steps done in the algorithm (if saved).
\item As well as other stuff for internal usage
}
}
\references{
Spiegel, Elmar, Thomas Kneib and Fabian Otto-Sobotka. Generalized additive models with flexible response functions. Statistics and Computing (2017). https://doi.org/10.1007/s11222-017-9799-6
}
\author{
Elmar Spiegel
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{predict.flexgam}}, \code{\link{plot.flexgam}}, \code{\link{deviance.flexgam}}, \code{\link{match_flexgam_control}}


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

plot(model_2, type = "response")
plot(model_2, type = "covariate")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
