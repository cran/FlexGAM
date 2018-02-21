\name{match_control}
\alias{match_flexgam_control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to check the control parameters
}
\description{
Controls if the values of the control parameters are correct. If any control parameters are missing those are set to default. For internal use only.
}
\usage{
match_flexgam_control(control = NULL, formula = formula, data=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{control}{
List of control parameters or NULL.
}
  \item{formula}{
Formula of the model to calculate the initial smoothing parameters.
}
  \item{data}{
Data of the model to calculate the initial smoothing parameters.
}
}
\details{
Controls if the values of the control parameters are correct. If control parameters are missing they are set to default:


\itemize{
    \item \code{"max_iter_in"}          =  \code{100}         Maximal number of inner iterations. 
    \item \code{"max_iter_out"}         =  \code{100}        Maximal number of outer iterations.
    \item \code{"delta_in"}             =  \code{1e-06}      Convergence of inner iterations. 
    \item \code{"delta_out"}            =  \code{1e-06}       Convergence of outer iterations. 
    \item \code{"min_mu_k"}             =  \code{1e-06}       Minimal value  of the fitted value. Also used to generate the upper limit for binomial data. Used to truncate the fitted values for numeric stability. (Occurrence can be read in the details). 
    \item \code{"min_Psi_d"}            =  \code{1e-06}      Minimal value of the derivative of the outer function. Used to truncate the derivatives for numeric stability. (Occurrence can be read in the details). 
    \item \code{"min_increase"}         =  \code{1e-04}       Minimal increase of the outer function.  
    %\item \code{"use_penalty"}          =  \code{TRUE}        Should the penalized deviance be used to check for non decreasing deviance? 
    %\item \code{"deriv_type"}           =  \code{"theoretic"} \code{theoretic}al or \code{numeric}al derivatives of the outer function? 
    %\item \code{"step_halving"}         =  \code{TRUE}        Should step-halving be applied 
    \item \code{"delta_halving"}        =  \code{1e-06}       Minimal difference at step-halving. 
    \item \code{"min_iter_halving_in"}  =  \code{1}           From which inner iteration should step halving be possible? 
    \item \code{"min_iter_halving_out"}  =  \code{2}          From which outer iteration the deviance stopping criterion should be applied? The minimum value is 2, to get the algorithm always starting.
    \item \code{"opt_function"}         =  \code{"optim"}    Which optimization function should be used to optimize the smoothing parameters? (\code{nlminb} or \code{optim(Nelder-Mead)}) 
    %\item \code{"start_logit"}          =  \code{TRUE}        Should the smoothing parameters of the standard \code{mgcv::gam} be used as initial values? 
    %\item \code{"start_grid_lambda"}    =  \code{TRUE}        Should a grid search be applied to get initial values for the smoothing parameter of the outer function? 
    \item \code{"initial_sm"}            =  \code{TRUE}        Should the smoothing parameters of the standard \code{mgcv::gam} be used as initial values for the covariates smoothing parameters and a grid search be applied to get initial values for the smoothing parameter of the outer function?
    \item \code{"fix_smooth"}           =  \code{FALSE}       Should the initial smoothing parameters (\code{sm_par_vec}) be used without optimization?  
    \item \code{"sm_par_vec"}           = \code{c("lambda"=1,"s(x1)"=\ldots)}    Initial smoothing parameters. Vector must start with \code{"lambda"} for the response function. The names of the covariate effects must fit to the \code{mgcv} output of the specified formula.  There is no need to specify the initial parameters, if \code{initial_sm = TRUE} and \code{fix_smooth = FALSE}.
    \item \code{"sp_range"}             =  \code{c(1e-8, 1e15)}   Range of all smoothing parameters. 
    %\item \code{"max_iter_out_opt"}   =  \code{25}          Maximal number of outer iterations during the optimization of the smoothing parameters.
    %\item \code{"max_iter_in_opt"}    =  \code{25}          Maximal number of inner iterations during the optimization of the smoothing parameters.
    \item \code{"reltol_opt"}         =  \code{1e-06}       Relative tolerance for optimizing the smoothing parameters. 
    %\item \code{"less_warnings"}        =  \code{TRUE}        Supress warnings during optimization of the smoothing parameters. 
    \item \code{"quietly"}              =  \code{FALSE}       Should the algorithm print steps of optimizing the smoothing parameters and iteration procedure for the final model? 
    \item \code{"save_step_response"}   =  \code{FALSE}       Should the steps of the algorithm be saved for convergences checks? 
    \item \code{"initial_model"}        =  \code{c("with_intercept","no_intercept")}       Whether the initial model should be estimated with or without intercept.
}


}

\value{
List of control parameters to fit the \code{flexgam} model.
}
\references{
Spiegel, Elmar, Thomas Kneib and Fabian Otto-Sobotka. Generalized additive models with flexible response functions. Statistics and Computing (2017). https://doi.org/10.1007/s11222-017-9799-6
}
\author{
Elmar Spiegel
}
\note{
The function is designed for internal usage.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{flexgam}}}
\examples{
# Only for internal usage.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
