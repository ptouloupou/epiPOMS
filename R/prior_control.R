#' Set hyperparameters for the prior distributions for epiPOMS_mcmc MCMC
#' algorithm.
#'
#' @description
#' Sets the hyperparameters for the prior distributions to be used in the
#' \code{\link{epiPOMS_mcmc}} Markov chain Monte Carlo (MCMC) algorithm to
#' perform Bayesian inference in partially observed multi-strain epidemic
#' models.
#'
#' @param
#' ntypes An integer greater than 1 providing the number of different strains in
#' the study, \eqn{n_g}.
#' @param
#' alphaprior A positive real number corresponding to the rate of the Exponential
#' prior distribution of the \eqn{\alpha} parameters, i.e.
#' the strain-specific external colonisation rates. Defaults to 1. See details,
#' below, for more information.
#' @param
#' betaprior A positive real number corresponding to the rate of the Exponential
#' prior distribution of the \eqn{\beta} parameters, i.e.
#' the strain-specific within-group colonisation rates. Defaults to 1. See
#' details, below, for more information.
#' @param
#' muprior A positive real number corresponding to the rate of the Exponential
#' prior distribution of the \eqn{\mu} parameters, i.e.
#' the strain-specific clearance rates. Defaults to 1. See details,
#' below, for more information.
#' @param
#' deltaprior A positive real number corresponding to the rate of the Exponential
#' prior distribution of the \eqn{\delta} parameter, i.e.
#' the relative colonisation rate in a carrier versus non-carrier individual.
#' Defaults to \eqn{\ln(2)}{ln(2)}. See details, below, for more information.
#' @param
#' gammaprior A positive real number corresponding to the rate of the Exponential
#' prior distribution of the \eqn{\gamma} parameter, i.e. the relative
#' colonisation rate in smaller versus bigger groups in terms of area (in square
#' meters). If there is no difference between groups, set \code{gammaprior =
#' NULL}. Defaults to \eqn{\ln(2)}{ln(2)}. See details, below, for more
#' information.
#' @param
#' nuprior A vector of positive real numbers providing the shape parameters of
#' the joint Dirichlet prior distribution of the \eqn{\nu} parameters, i.e.
#' the probabilities of carriage at the beginning of the study. This is a vector
#' of \eqn{n_g + 1} values. The first element corresponds to the non-carriage
#' state. The second and subsequent elements correspond to the carriage of one of
#' the \eqn{n_g} strains.  The default value for all \code{nuprior}
#' elements is 0.5. See details, below, for more information.
#' @param
#' theta12prior A vector of two positive real numbers providing the shape
#' parameters of the Beta prior distribution of the \eqn{\theta_1} and
#' \eqn{\theta_2} parameters, i.e. the test sensitivities. By default,
#' \code{theta12prior = c(0.5, 0.5)}. See details, below, for more
#' information.
#' @param
#' thetaCSprior A vector of three positive real numbers providing the shape
#' parameters of the joint Dirichlet prior distribution of the \eqn{\theta_C} and
#' \eqn{\theta_S} parameters, i.e. given that a test is found positive
#' \eqn{\theta_C} denotes the probability of correctly identifying a
#' common strain \eqn{(1, 2, \ldots, n_g - 1)} and \eqn{\theta_S} is the
#' probability of misclassifying a common strain with a different common strain.
#' The default values for the \code{thetaCSprior} elements are (4.5, 0.25, 0.25).
#' See details, below, for more information.
#' @param
#' thetaPprior A vector of two positive real numbers providing the shape
#' parameters of the Beta prior distribution of the \eqn{\theta_P} parameter,
#' i.e. the probability that a strain of pooled type \eqn{n_g} is
#' classified as a common strain given that a test is found positive.
#' By default, \code{thetaPprior = c(0.5, 0.5)}. See details, below, for more
#' information.
#'
#' @return
#' A list with components \code{alphaprior}, \code{betaprior}, \code{muprior},
#' \code{deltaprior}, \code{gammaprior}, \code{nuprior}, \code{theta12prior},
#' \code{thetaCSprior} and \code{thetaPprior}. This can be used as an argument
#' of functions \code{\link{epiPOMS_mcmc}} and \code{\link{mcmc_initpars}}.
#'
#' @details
#' Auxiliary function that can be used to set hyperparameter values for the
#' prior distributions for the MCMC algorithm used by \code{\link{epiPOMS_mcmc}}
#' for performing Bayesian inference in partially observed multi-strain epidemic
#' models.
#'
#' For the strain-specific external colonisation rates (\eqn{\alpha}), the
#' within-pen colonisation rates (\eqn{\beta}) and the clearance rates
#' (\eqn{\mu}), we assign univariate Exponential priors each with default rate
#' equal to 1. The priors for \eqn{\delta} and \eqn{\gamma} are also assumed to
#' be Exponential with default rate equal to \eqn{\ln(2)}{ln(2)}, reflecting
#' equal prior probabilities for these parameters to be less or more than
#' one. By default parameter \eqn{\delta} is included in the model. However,
#' if there is no difference between groups, we allow for the possibility to
#' exclude this parameter for the model by setting \code{gammaprior = NULL}.
#' For the misspecification parameters \eqn{\theta_1}, \eqn{\theta_2} and
#' \eqn{\theta_P} we assume Beta prior distributions, with default
#' of Beta(0.5, 0.5) which is the Jeffreys' prior
#' \insertCite{Jeffreys1961}{epiPOMS}. For the remaining observation
#' parameters we assume a joint Dirichlet prior distribution, with default the
#' \eqn{(\theta_C,} \eqn{\theta_S,} \eqn{1 - \theta_C - \theta_S)}  \eqn{\sim}{~}
#' Dirichlet(4.5, 0.25, 0.25) minimally informative prior
#' \insertCite{Kelly2011}{epiPOMS}. Finally, for the
#' probabilities of carriage at the beginning of the study (\eqn{\nu}) we use
#' a joint Dirichlet prior with a default the Jeffreys' Dirichlet
#' distribution with all \eqn{n_g + 1} parameters set to 0.5. For the full
#' details see \insertCite{Touloupou2020;textual}{epiPOMS}.
#'
#' @author
#' Panayiota Touloupou
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for performing inference in partially observed
#' multi-strain epidemic models and \code{\link{mcmc_initpars}} for producing
#' the initial values for the MCMC algorithm.
#'
#' @export
prior_control <- function(ntypes, alphaprior = 1, betaprior = 1, muprior = 1, deltaprior = log(2),
    gammaprior = log(2), nuprior = rep(0.5, ntypes + 1), theta12prior = c(0.5, 0.5), thetaCSprior = c(4.5,
        0.25, 0.25), thetaPprior = c(0.5, 0.5)) {
    
    # Error checks for input arguments
    check <- ntypes
    if (is.null(check)) {
        stop("prior_control: The ntypes has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("prior_control: The ntypes must be an integer greater than one.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("prior_control: The ntypes must be an integer greater than one.", call. = FALSE)
    } else if ((check <= 1) | (check != round(ntypes))) {
        stop("prior_control: The ntypes must be an integer greater than one.", call. = FALSE)
    }
    
    check <- alphaprior
    if (is.null(check)) {
        stop("prior_control: The alphaprior has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("prior_control: The alphaprior must be a positive real number.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("prior_control: The alphaprior must be a positive real number.", call. = FALSE)
    } else if ((check <= 0)) {
        stop("prior_control: The alphaprior must be a positive real number.", call. = FALSE)
    }
    
    check <- betaprior
    if (is.null(check)) {
        stop("prior_control: The betaprior has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("prior_control: The betaprior must be a positive real number.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("prior_control: The betaprior must be a positive real number.", call. = FALSE)
    } else if ((check <= 0)) {
        stop("prior_control: The betaprior must be a positive real number.", call. = FALSE)
    }
    
    check <- muprior
    if (is.null(check)) {
        stop("prior_control: The muprior has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("prior_control: The muprior must be a positive real number.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("prior_control: The muprior must be a positive real number.", call. = FALSE)
    } else if ((check <= 0)) {
        stop("prior_control: The muprior must be a positive real number.", call. = FALSE)
    }
    
    check <- deltaprior
    if (is.null(check)) {
        stop("prior_control: The deltaprior has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("prior_control: The deltaprior must be a positive real number.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("prior_control: The deltaprior must be a positive real number.", call. = FALSE)
    } else if ((check <= 0)) {
        stop("prior_control: The deltaprior must be a positive real number.", call. = FALSE)
    }
    
    check <- gammaprior
    if (!(is.null(check))) {
        if ((length(check) != 1) | is.list(check)) {
            stop("prior_control: The gammaprior must be a positive real number or NULL.", call. = FALSE)
        } else if (!is.numeric(check)) {
            stop("prior_control: The gammaprior must be a positive real number or NULL.", call. = FALSE)
        } else if ((check <= 0)) {
            stop("prior_control: The gammaprior must be a positive real number or NULL.", call. = FALSE)
        }
    }
    
    check <- nuprior
    if (is.null(check)) {
        stop("prior_control: The nuprior has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("prior_control: The nuprior must be a vector.", call. = FALSE)
    } else if (length(check) != (ntypes + 1)) {
        stop("prior_control: Length of nuprior is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("prior_control: The entries of nuprior must be positive real numbers.", call. = FALSE)
    } else if (any(check <= 0)) {
        stop("prior_control: The entries of nuprior must be positive real numbers.", call. = FALSE)
    }
    
    check <- theta12prior
    if (is.null(check)) {
        stop("prior_control: The theta12prior has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("prior_control: The theta12prior must be a vector.", call. = FALSE)
    } else if (length(check) != 2) {
        stop("prior_control: Length of theta12prior is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("prior_control: The entries of theta12prior must be positive real numbers.", call. = FALSE)
    } else if (any(check <= 0)) {
        stop("prior_control: The entries of theta12prior must be positive real numbers.", call. = FALSE)
    }
    
    check <- thetaCSprior
    if (is.null(check)) {
        stop("prior_control: The thetaCSprior has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("prior_control: The thetaCSprior must be a vector.", call. = FALSE)
    } else if (length(check) != 3) {
        stop("prior_control: Length of thetaCSprior is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("prior_control: The entries of thetaCSprior must be positive real numbers.", call. = FALSE)
    } else if (any(check <= 0)) {
        stop("prior_control: The entries of thetaCSprior must be positive real numbers.", call. = FALSE)
    }
    
    check <- thetaPprior
    if (is.null(check)) {
        stop("prior_control: The thetaPprior has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("prior_control: The thetaPprior must be a vector.", call. = FALSE)
    } else if (length(check) != 2) {
        stop("prior_control: Length of thetaPprior is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("prior_control: The entries of thetaPprior must be positive real numbers.", call. = FALSE)
    } else if (any(check <= 0)) {
        stop("prior_control: The entries of thetaPprior must be positive real numbers.", call. = FALSE)
    }
    # End of error checks for input arguments
    
    output <- list(alphaprior = alphaprior, betaprior = betaprior, muprior = muprior, deltaprior = deltaprior, 
        gammaprior = gammaprior, nuprior = nuprior, theta12prior = theta12prior, thetaCSprior = thetaCSprior, 
        thetaPprior = thetaPprior)
    
    return(output)
}
