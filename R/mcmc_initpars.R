#' Set starting values for epiPOMS_mcmc MCMC algorithm.
#'
#' @description
#' Sets starting values used to initialise the Markov chain Monte Carlo (MCMC)
#' algorithm by \code{\link{epiPOMS_mcmc}} function.
#'
#' @param
#' inittype A character string determining which of the two options to compute the
#' initial values for the model parameters should be used. Options are
#' "prior" (the default) and "manual". Setting \code{inittype = } \code{"prior"}
#' initial values are drawn from their prior distribution; \code{inittype = }
#' \code{"manual"} initial values are provided by the user. See details, below,
#' for more information.
#' @param
#' priorcontrol A list containing the hyperparameters for the prior
#' distributions for the MCMC algorithm, obtained using function
#' \code{\link{prior_control}}. Only used if \code{inittype = } \code{"prior"}.
#' @param
#' initpars A list containing the starting values of the model parameters to be
#' used to initialise the MCMC algorithm, if \code{inittype = } \code{"manual"}.
#' The names of the components must be "alpha", "beta", "mu", "delta", "gamma",
#' "nu" and "theta", where:
#' \describe{
#' \item{alpha}{is a vector of non-negative real numbers, one for each strain,
#' containing the initial values of the \eqn{\alpha} parameters, i.e. the
#' strain-specific external colonisation rates. This is a vector of \eqn{n_g}
#' values, where \eqn{n_g} is the number of different strains in the study.
#' The first element corresponds to the first strain, the second element to
#' the second strain and so on.}
#' \item{beta}{is a vector of non-negative real numbers, one for each strain,
#' containing the initial values of the \eqn{\beta} parameters, i.e. the
#' strain-specific within-group colonisation rates. This is a vector of
#' \eqn{n_g} values. The first element corresponds to the first strain, the
#' second element to the second strain and so on.}
#' \item{mu}{is a vector of non-negative real numbers, one for each strain,
#' containing the initial values of the \eqn{\mu} parameters, i.e. the
#' strain-specific clearance rates. This is a vector of \eqn{n_g} values.
#' The first element corresponds to the first strain, the second element to
#' the second strain and so on.}
#' \item{delta}{is a non-negative real number corresponding to the initial value
#' of the \eqn{\delta} parameter, i.e. the relative colonisation rate in a
#' carrier versus non-carrier individual.}
#' \item{gamma}{is a non-negative real number corresponding to the initial value
#' of the \eqn{\gamma} parameter, i.e. the relative colonisation rate in smaller
#' versus bigger groups in terms of area (in square meters). If
#' \code{smallgroups = NULL}, i.e there is no difference between groups, set
#' \code{gamma = NULL}.}
#' \item{nu}{is a vector of numbers between 0 and 1, containing the initial values
#' of the \eqn{\nu} parameters, i.e. the probabilities of carriage at the
#' beginning of the study. This is a vector of \eqn{(n_g + 1)} values. The first
#' element corresponds to the non-carriage state. The second and subsequent
#' elements correspond to the carriage of one of the \eqn{n_g} strains. These
#' probabilities should sum up to one.}
#' \item{theta}{is a vector of 5 values, containing the initial values for the
#' following parameters in order: \eqn{\theta_1}, \eqn{\theta_2}, \eqn{\theta_C},
#' \eqn{\theta_S} and \eqn{\theta_P}, where \eqn{\theta_1} and \eqn{\theta_2}
#' denote the test sensitivities. Given that a test is positive, \eqn{\theta_C}
#' denotes the probability of correctly identifying a common strain,
#' \eqn{\theta_S} is the probability of misclassifying a common strain with a
#' different common strain and \eqn{\theta_P} the probability that a strain of
#' pooled type is classified as a common strain. All parameters must be between
#' 0 and 1 and \eqn{\theta_C + \theta_S \leq 1}{\theta_C + \theta_S \le 1}.}
#' }
#' @param
#' epidata An object of class \sQuote{epiPOMSdata}, produced from the
#' \code{\link{as_epiPOMSdata}} function.
#' @param
#' smallgroups A vector of positive integers, with length equal to the number
#' of groups that are smaller in terms of area (in square meters), containing
#' the indices of the small groups. These values must be between 1 and the total
#' number of groups in the study. Set \code{smallgroups = NULL} if there is no
#' difference between groups.
#'
#' @return
#' A list is returned with the following components:
#' \describe{
#' \item{init_alpha}{A vector of initial values of the \eqn{\alpha}
#' parameters, i.e. the strain-specific external colonisation rates.
#' This is a vector of \eqn{n_g} values, where \eqn{n_g} is the number of
#' strains in the study. The first element corresponds to the first strain,
#' the second element to the second strain and so on.}
#' \item{init_beta}{A vector of initial values of the \eqn{\beta}
#' parameters, i.e. the strain-specific within-group colonisation rates.
#' This is a vector of \eqn{n_g} values. The first element corresponds to the
#' first strain, the second element to the second strain and so on.}
#' \item{init_mu}{A vector of initial values of the \eqn{\mu}
#' parameters, i.e. the strain-specific clearance rates.
#' This is a vector of \eqn{n_g} values. The first element corresponds to the
#' first strain, the second element to the second strain and so on.}
#' \item{init_delta}{Initial value of the \eqn{\delta} parameter, i.e.
#' the relative colonisation rate in a carrier versus non-carrier individual.}
#' \item{init_gamma}{Initial value of the \eqn{\gamma} parameter, i.e.
#' the relative colonisation rate in smaller versus bigger groups in terms of
#' area (in square meters), if \code{smallgroups} \eqn{\neq}{≠} \code{NULL}.
#' Otherwise, \code{init_gamma = NULL}.}
#' \item{init_nu}{A vector of initial values of the \eqn{\nu} parameters,
#' i.e. the probabilities of carriage at the beginning of the study. This is a
#' vector of \eqn{(n_g + 1)} values. The first element corresponds to the
#' non-carriage state. The second and subsequent elements correspond to the
#' carriage of one of the \eqn{n_g} strains.}
#' \item{init_theta}{A vector of initial values of the observation
#' parameters, in order, \eqn{\theta_1}, \eqn{\theta_2}, \eqn{\theta_C},
#' \eqn{\theta_S}, \eqn{\theta_P}, where \eqn{\theta_1} and \eqn{\theta_2}
#' denote the test sensitivities. Given that a test is positive, \eqn{\theta_C}
#' denotes the probability of correctly identifying a common strain,
#' \eqn{\theta_S} is the probability of misclassifying a common strain with a
#' different common strain and \eqn{\theta_P} the probability that a strain of
#' pooled type is classified as a common strain.}
#' \item{init_hidprocess}{A list containing \eqn{P} matrices (one for each
#' group), with the initial hidden carriage process, where \eqn{P} is the total
#' number of groups in the study. For more information, see
#' \code{\link{hidstate_sim}}.}
#' }
#'
#' @details
#' Auxiliary function that can be used to generate initial values for the MCMC
#' algorithm used by \code{\link{epiPOMS_mcmc}} for performing Bayesian inference
#' in partially observed multi-strain epidemic models. The output of this function
#' can be used as an argument of \code{\link{mcmc_control}}.
#'
#' The argument \code{inittype} has two options. When \code{inittype = }
#' \code{"prior"} the initial values are generated by sampling each model
#' parameter from its prior distribution. In this case, the argument
#' \code{priorcontrol} has to be specified. When  \code{inittype = }
#' \code{"manual"} the initial values are provided by the user, through
#' argument \code{initpars}, which has to be specified. The default is to sample
#' the initial values directly from the prior distribution of the model
#' parameters. However, in some cases (when using very vague or even improper
#' flat priors), this is not a good idea since it could generate disparate
#' initial values that could make the chain take an excessive time to achieve
#' convergence. Therefore, as a general guideline, it is a safe procedure to give
#' the initial values of the parameters, at least for those having vague or
#' improper distributions.
#'
#' @author
#' Panayiota Touloupou
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for performing inference in partially observed
#' multi-strain epidemic models, \code{\link{prior_control}} for specifying the
#' hyperparameters of the prior distributions and \code{\link{mcmc_control}}
#' for setting the parameter values that control the MCMC algorithm.
#'
#' @export
mcmc_initpars <- function(inittype = "prior", priorcontrol, initpars, epidata, smallgroups) {
    
    # Error checks for input arguments
    
    check <- inittype
    if (is.null(check)) {
        stop("mcmc_initpars: The inittype has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("mcmc_initpars: The inittype must be a character string, either  'prior' or 'manual'.", 
            call. = FALSE)
    } else if (!(check %in% c("prior", "manual"))) {
        stop("mcmc_initpars: The inittype must be a character string, either  'prior' or 'manual'.", 
            call. = FALSE)
    }
    
    if (inittype == "prior") {
        check <- priorcontrol
        checknames <- c("alphaprior", "betaprior", "muprior", "deltaprior", "gammaprior", "nuprior", 
            "theta12prior", "thetaCSprior", "thetaPprior")
        if (is.null(check)) {
            stop("mcmc_initpars: The priorcontrol has to be specified.", call. = FALSE)
        } else if (!is.list(check)) {
            stop("mcmc_initpars: The priorcontrol must be a list obtained by function prior_control.", 
                call. = FALSE)
        } else if (length(check) != 9) {
            stop("mcmc_initpars: The priorcontrol must be a list obtained by function prior_control.", 
                call. = FALSE)
        } else if (any(names(check) != checknames)) {
            stop("mcmc_initpars: The priorcontrol must be a list obtained by function prior_control.", 
                call. = FALSE)
        }
        
        if ((is.null(smallgroups)) != (is.null(check$gammaprior))) {
            stop("mcmc_initpars: The smallgroups must be NULL if gammaprior at priorcontrol is NULL, and via versa.", 
                call. = FALSE)
        }
        
    } else {
        
        infoepidata <- info_epiPOMSdata(epidata)
        
        check <- initpars$alpha
        if (is.null(check)) {
            stop("mcmc_initpars: The alpha at initpars has to be specified.", call. = FALSE)
        } else if ((!is.vector(check)) | is.list(check)) {
            stop("mcmc_initpars: The alpha at initpars must be a vector.", call. = FALSE)
        } else if (length(check) != infoepidata$ntypes) {
            stop("mcmc_initpars: Length of alpha at initpars is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("mcmc_initpars: The entries of alpha at initpars must be non-negative reals.", 
                call. = FALSE)
        } else if (any(check < 0)) {
            stop("mcmc_initpars: The entries of alpha at initpars must be non-negative reals.", 
                call. = FALSE)
        }
        
        check <- initpars$beta
        if (is.null(check)) {
            stop("mcmc_initpars: The beta at initpars has to be specified.", call. = FALSE)
        } else if ((!is.vector(check)) | is.list(check)) {
            stop("mcmc_initpars: The beta at initpars must be a vector.", call. = FALSE)
        } else if (length(check) != infoepidata$ntypes) {
            stop("mcmc_initpars: Length of beta at initpars is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("mcmc_initpars: The entries of beta at initpars must be non-negative reals.", 
                call. = FALSE)
        } else if (any(check < 0)) {
            stop("mcmc_initpars: The entries of beta at initpars must be non-negative reals.", 
                call. = FALSE)
        }
        
        check <- initpars$mu
        if (is.null(check)) {
            stop("mcmc_initpars: The mu at initpars has to be specified.", call. = FALSE)
        } else if ((!is.vector(check)) | is.list(check)) {
            stop("mcmc_initpars: The mu at initpars must be a vector.", call. = FALSE)
        } else if (length(check) != infoepidata$ntypes) {
            stop("mcmc_initpars: Length of mu at initpars is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("mcmc_initpars: The entries of mu at initpars must be non-negative reals.", call. = FALSE)
        } else if (any(check < 0)) {
            stop("mcmc_initpars: The entries of mu at initpars must be non-negative reals.", call. = FALSE)
        }
        
        check <- initpars$delta
        if (is.null(check)) {
            stop("mcmc_initpars: The delta at initpars has to be specified.", call. = FALSE)
        } else if ((length(check) != 1) | is.list(check)) {
            stop("mcmc_initpars: The delta at initpars must be a non-negative real.", call. = FALSE)
        } else if (!is.numeric(check)) {
            stop("mcmc_initpars: The delta at initpars must be a non-negative real.", call. = FALSE)
        } else if ((check < 0)) {
            stop("mcmc_initpars: The delta at initpars must be a non-negative real.", call. = FALSE)
        }
        
        check <- initpars$gamma
        if (!is.null(check)) {
            if ((length(check) != 1) | is.list(check)) {
                stop("mcmc_initpars: The gamma at initpars must be a non-negative real or NULL.", 
                  call. = FALSE)
            } else if (!is.numeric(check)) {
                stop("mcmc_initpars: The gamma at initpars must be a non-negative real or NULL.", 
                  call. = FALSE)
            } else if ((check < 0)) {
                stop("mcmc_initpars: The gamma at initpars must be a non-negative real or NULL.", 
                  call. = FALSE)
            }
        }
        
        if ((is.null(smallgroups)) != (is.null(check))) {
            stop("mcmc_initpars: The smallgroups must be NULL if gamma at initpars is NULL, and via versa.", 
                call. = FALSE)
        }
        
        
        check <- initpars$nu
        if (is.null(check)) {
            stop("mcmc_initpars: The nu at initpars has to be specified.", call. = FALSE)
        } else if ((!is.vector(check)) | is.list(check)) {
            stop("mcmc_initpars: The nu at initpars must be a vector.", call. = FALSE)
        } else if (length(check) != (infoepidata$ntypes + 1)) {
            stop("mcmc_initpars: Length of nu at initpars is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("mcmc_initpars: The entries of nu at initpars must be between zero and one.", 
                call. = FALSE)
        } else if (any(check < 0) | any(check > 1)) {
            stop("mcmc_initpars: The entries of nu at initpars must be between zero and one.", 
                call. = FALSE)
        } else if (sum(check) > 1) {
            stop("mcmc_initpars:  The entries of nu at initpars must sum to one.", call. = FALSE)
        }
        
        check <- initpars$theta
        if (is.null(check)) {
            stop("mcmc_initpars: The theta at initpars has to be specified.", call. = FALSE)
        } else if ((!is.vector(check)) | is.list(check)) {
            stop("mcmc_initpars: The theta at initpars must be a vector.", call. = FALSE)
        } else if (length(check) != 5) {
            stop("mcmc_initpars: Length of theta at initpars is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("mcmc_initpars: The entries of theta at initpars must be between zero and one.", 
                call. = FALSE)
        } else if (any(check < 0) | any(check > 1)) {
            stop("mcmc_initpars: The entries of theta at initpars must be between zero and one.", 
                call. = FALSE)
        } else if (sum(check[3:4]) > 1) {
            stop("mcmc_initpars: The sum of third and forth elements of theta at initpars must be less than one.", 
                call. = FALSE)
        }
        
    }
    if (!is(epidata, "epiPOMSdata")) {
        stop("mcmc_initpars: The epidata must be a class of 'epiPOMSdata'.", call. = FALSE)
    }
    
    infoepidata <- info_epiPOMSdata(epidata)
    
    check <- smallgroups
    if (!is.null(check)) {
        if ((!is.vector(check)) | is.list(check)) {
            stop("mcmc_initpars: The smallgroups must be a vector or NULL.", call. = FALSE)
        } else if (length(check) > infoepidata$ngroups) {
            stop("mcmc_initpars: Length of smallgroups is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("mcmc_initpars: The entries of smallgroups must be integers between 1 and number of groups in the study or NULL.", 
                call. = FALSE)
        } else if (any(check < 1) | any(check > infoepidata$ngroups) | any(check != round(check))) {
            stop("mcmc_initpars: The entries of smallgroups must be integers between 1 and number of groups in the study or NULL.", 
                call. = FALSE)
        }
    }
    
    
    # End of error checks for input arguments
    
    infoepidata <- info_epiPOMSdata(epidata)
    
    if (inittype == "prior") {
        
        init_alpha <- rexp(n = infoepidata$ntypes, rate = priorcontrol$alphaprior)
        init_beta <- rexp(n = infoepidata$ntypes, rate = priorcontrol$betaprior)
        init_mu <- rexp(n = infoepidata$ntypes, rate = priorcontrol$muprior)
        init_delta <- rexp(n = 1, rate = priorcontrol$deltaprior)
        if (!is.null(priorcontrol$gammaprior)) {
            init_gamma <- rexp(n = 1, rate = priorcontrol$gammaprior)
        } else {
            init_gamma <- NULL
        }
        init_nu <- as.vector(MCMCpack::rdirichlet(n = 1, alpha = priorcontrol$nuprior))
        if (sum(init_nu) != 1) {
            init_nu[1] <- 1 - sum(init_nu[-1])
        }
        
        init_theta <- rep(0, 5)
        init_theta[1:2] <- rbeta(n = 2, shape1 = priorcontrol$theta12prior[1], shape2 = priorcontrol$theta12prior[2])
        init_theta[3:4] <- MCMCpack::rdirichlet(n = 1, alpha = priorcontrol$thetaCSprior)[1:2]
        init_theta[5] <- rbeta(n = 1, shape1 = priorcontrol$thetaPprior[1], shape2 = priorcontrol$thetaPprior[2])
        
    } else {
        
        init_alpha <- initpars$alpha
        init_beta <- initpars$beta
        init_mu <- initpars$mu
        init_delta <- initpars$delta
        init_gamma <- initpars$gamma
        init_nu <- initpars$nu
        init_theta <- initpars$theta
        
    }
    
    init_hidpars <- list(alpha = init_alpha, beta = init_beta, mu = init_mu, delta = init_delta, 
        gamma = init_gamma, nu = init_nu)
    
    init_hidprocess <- hidstate_sim(ntypes = infoepidata$ntypes, ngroups = infoepidata$ngroups, 
        ninds = infoepidata$ninds, tmax = infoepidata$tmax, hidpars = init_hidpars, smallgroups = smallgroups, 
        indNA = infoepidata$indNA, tNA = unlist(lapply(infoepidata$tmaxNA, min)))
    
    classmatrices <- classprob(ntypes = infoepidata$ntypes, obserpars = init_theta)
    
    for (i in 1:infoepidata$ngroups) {
        gamma_new <- ifelse(i %in% smallgroups, init_gamma, 1)
        init_hidpars$gamma <- gamma_new
        
        for (j in 1:infoepidata$ninds[i]) {
            init_hidprocess[[i]][j, ] <- iFFBS_upd(ntypes = infoepidata$ntypes, nindsgr = infoepidata$ninds[i], 
                ind = j, tmaxgr = infoepidata$tmax[[i]], tmaxind = infoepidata$tmaxNA[[i]][j], 
                curhidprocgr = init_hidprocess[[i]], obsertest1ind = infoepidata$obserdata$test1[[i]][j, 
                  ], obsertest2ind = infoepidata$obserdata$test2[[i]][j, ], curalpha = init_hidpars$alpha, 
                curbeta = init_hidpars$beta, curmu = init_hidpars$mu, curdelta = init_hidpars$delta, 
                curgamma = init_hidpars$gamma, curnu = init_hidpars$nu, curclassmat = classmatrices)
        }
    }
    
    
    output <- list(init_alpha = init_alpha, init_beta = init_beta, init_mu = init_mu, init_delta = init_delta, 
        init_gamma = init_gamma, init_nu = init_nu, init_theta = init_theta, init_hidprocess = init_hidprocess)
    
    return(output)
}



