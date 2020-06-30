#' Set control parameters for epiPOMS_mcmc MCMC algorithm.
#'
#' @description
#' Sets parameter values that control the Markov chain Monte Carlo (MCMC)
#' algorithm used by \code{\link{epiPOMS_mcmc}} function to perform
#' Bayesian inference in partially observed multi-strain epidemic models.
#'
#' @param
#' nsamp A positive integer giving the number of iterations to run the MCMC
#' algorithm. Defaults to 35000.
#' @param
#' burnin A non-negative integer giving the number of burn-in iterations to be
#' run by the MCMC algorithm, with \code{burnin} \eqn{<} \code{nsamp}. Defaults
#' to 10000.
#' @param
#' thinning A positive integer corresponding to thinning parameter; 1 in
#' \code{thinning} iterations will be kept, after the \code{burnin} phase,
#' the rest will be discarded. Requirements: \code{thinning} \eqn{\leq}{\le}
#' (\code{nsamp} - \code{burnin}). Defaults to 5.
#' @param
#' seed An integer used as the seed for the random number generator (RNG) at
#' the start of the MCMC algorithm. By default, \code{seed = NULL}, which
#' does not alter the RNG state.
#' @param
#' init_pars A list containing the initial values for the MCMC algorithm,
#' obtained using function \code{\link{mcmc_initpars}}.
#' @param
#' init_stepsize A positive real value specifying the initial leapfrog
#' stepsize in Hamiltonian Monte Carlo (HMC) method. Defaults to 0.05.
#' See details, below, for more information.
#' @param
#' max_nsteps A positive integer corresponding to the maximum number of leapfrog
#' steps in HMC method. Defaults to 20. See details, below, for more
#' information.
#' @param
#' ntypes An integer greater than 1 providing the number of different strains in
#' the study, \eqn{n_g}.
#' @param
#' mass_diag A vector of positive real numbers providing the diagonal of the mass
#' matrix, with length equal to the number of parameters to be updated jointly
#' using the HMC algorithm. The length of the vector should be
#' \eqn{(3 \times n_g + 2)}{(3xn_g + 2)} if parameter \eqn{\gamma} is included
#' in the model, otherwise it should be \eqn{(3 \times n_g + 1)}{(3xn_g + 1)}.
#' See \code{\link{HMC_upd}} for more details. By default, all elements of
#' \code{mass_diag} are equal to one.
#'
#' @return
#' A list with components \code{nsamp}, \code{burnin}, \code{thinning},
#' \code{seed}, \code{init_pars}, \code{init_stepsize}, \code{max_nsteps}
#' and \code{mass_diag}. This can be used as an argument of function
#' \code{\link{epiPOMS_mcmc}}.
#'
#' @details
#' Auxiliary function that can be used to set parameter values that control
#' the MCMC algorithm used by \code{\link{epiPOMS_mcmc}} for performing Bayesian
#' inference in multi-strain epidemic models. This function is only used in
#' conjunction with the \code{\link{epiPOMS_mcmc}} function.
#'
#' Sampling from the posterior distribution is done by using an MCMC algorithm
#' that employs both Gibbs and HMC techniques (see
#' \insertCite{Touloupou2020;textual}{epiPOMS} for more details). The MCMC
#' scheme starts by setting initial parameter values or drawing them from
#' their prior distribution, \code{init_pars}, as returned by function
#' \code{\link{mcmc_initpars}}. The MCMC is then run for a total number of
#' \code{nsamp} iterations but the output is only recorded after the
#' \code{burnin}, and only 1 in every \code{thinning} iterations will be kept,
#' so that the posterior sample size for each of the model parameters is
#' \code{round}((\code{nsamp} - \code{burnin})/\code{thinning}).
#'
#' In the sampling process of HMC, leapfrog integration is used to
#' approximately solve Hamilton's equations, and the integration is required
#' to set the number of discrete time steps and the integration stepsize. These
#' are tuning parameters for the HMC algorithm. The approach adopted in our
#' study is to randomly select the number of leapfrog steps at each iteration of
#' HMC, uniformly from (0, \code{max_nsteps}) and adapt the stepsize during
#' burnin to obtain an acceptance rate of roughly 65\% as suggested by
#' \insertCite{Neal2011;textual}{epiPOMS}, with a starting value of
#' \code{init_stepsize}. The mass matrix should be symmetric and
#' positive-definite, which is typically diagonal, and is often a scalar
#' multiple of the identity matrix. Here, we assume that the mass matrix is
#' diagonal, with diagonal elements being the elements of \code{mass_diag}.
#'
#' The default values work well in the MCMC algorithm for the analysis of the
#' real dataset concerning the transmission dynamics of \emph{Escherichia coli}
#' O157:H7 in cattle, see \code{\link{Ecoliepidata}}. For the analysis of other
#' infectious diseases or datasets, these values may need to be adjusted by the
#' user to achieve good mixing properties.
#'
#' @author
#' Panayiota Touloupou
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for performing inference in partially observed
#' multi-strain epidemic  models and \code{\link{mcmc_initpars}} for
#' producing the initial values for the MCMC algorithm.
#'
#' @examples
#' # Simulate partially observed multi-strain epidemic data. Individuals in
#' # 15 different groups of size varied from 4 to 9, were tested at the same
#' # time points, every 3 days during a period of 100 days. Five strains
#' # were identified.
#' set.seed(1)
#' ntypes <- 5
#' ngroups <- 15
#' ninds <- sample(4:9, size = ngroups, replace = TRUE)
#' tmax <- rep(100, ngroups)
#' hidpars <- list(alpha = c(sample(c(0.0015, 0.001), ntypes - 1, replace = TRUE),
#'      0.002), beta = sample(c(0.005, 0.01), ntypes, replace = TRUE), mu =
#'      sample(c(0.1, 0.15), ntypes, replace = TRUE), delta = 0.5, gamma = 2,
#'      nu = c(0.9, rep((1-0.9)/(ntypes+1), (ntypes-1)), 2*(1-0.9)/(ntypes+1)))
#' obserpars <- c(0.8, 0.5, 0.8, 0.025, 0.05)
#' smallgroups <- sample(1:ngroups, size = 10)
#' obsertimes <- lapply(tmax, function(x) seq(1, x, 3))
#' indNA <- c(rep(0, ngroups-1), 1)
#' tNA <- c(rep(100, ngroups-1), 70)
#' npostyped <- 20
#' epidata <- obserdata_sim(ntypes = ntypes, ngroups = ngroups, ninds = ninds,
#'      tmax = tmax, hidpars = hidpars, obserpars = obserpars, smallgroups =
#'      smallgroups, obsertimes = obsertimes, indNA = indNA, tNA = tNA,
#'      npostyped = npostyped)
#'
#' # Modify and extract information
#' infoepidata <- info_epiPOMSdata(epidata)
#'
#' # Set hyperparameters for the prior distributions
#' priorcontrol <- prior_control(ntypes = ntypes)
#'
#' # Set initial values for MCMC algorithm
#' initpars <- list(alpha = rep(0.01, infoepidata$ntypes), beta = rep(0.01,
#'      infoepidata$ntypes), mu = rep(0.01, infoepidata$ntypes), delta = 1,
#'      gamma = 1, nu =  c(0.9, rep(0.1/(infoepidata$ntypes),
#'      infoepidata$ntypes)), theta = c(0.5, 0.5, 0.9, 0.05, 0.5))
#' mcmcinitpars <- mcmc_initpars(inittype = "manual", initpars = initpars,
#'      epidata = epidata, smallgroups = smallgroups)
#'
#' # Set control parameters for MCMC algorithm
#' mcmccontrol <- mcmc_control(nsamp = 50, burnin = 10, init_pars =
#'      mcmcinitpars, ntypes = ntypes)
#'
#' # Run MCMC algorithm on the simulated epidemic
#' # Note: Not enough iterations for any real inference
#' # Note: This will take a few seconds to run
#' tic = Sys.time()
#' examplemcmc <- epiPOMS_mcmc(infoepidata = infoepidata, mcmccontrol =
#'      mcmccontrol, priorcontrol = priorcontrol, smallgroups = smallgroups,
#'      parallel = FALSE)
#' toc = Sys.time() - tic
#' toc
#'
#' # Note: Not enough iterations for assessing convergence or mixing of chains
#' plot(x = examplemcmc, infoepidata = infoepidata)
#'
#' @export
mcmc_control <- function(nsamp = 35000, burnin = 10000, thinning = 5, seed = NULL, init_pars, init_stepsize = 0.05, 
    max_nsteps = 20, ntypes, mass_diag = if (is.null(init_pars$init_gamma)) rep(1, 3 * ntypes + 
        1) else rep(1, 3 * ntypes + 2)) {
    
    # Error checks for input arguments
    check <- nsamp
    if (is.null(check)) {
        stop("mcmc_control: The nsamp has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("mcmc_control: The nsamp must be a positive integer.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("mcmc_control: The nsamp must be a positive integer.", call. = FALSE)
    } else if ((check <= 0) | (check != round(check))) {
        stop("mcmc_control: The nsamp must be a positive integer.", call. = FALSE)
    }
    
    check <- burnin
    if (is.null(check)) {
        stop("mcmc_control: The burnin has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("mcmc_control: The burnin must be a non-negative integer less than nsamp.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("mcmc_control: The burnin must be a non-negative integer less than nsamp.", call. = FALSE)
    } else if ((check < 0) | (check >= nsamp) | (check != round(check))) {
        stop("mcmc_control: The burnin must be a non-negative integer less than nsamp.", call. = FALSE)
    }
    
    check <- thinning
    if (is.null(check)) {
        stop("mcmc_control: The thinning has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("mcmc_control: The thinning must be a positive integer equal or less than (nsamp - burnin).", 
            call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("mcmc_control: The thinning must be a positive integer equal or less than (nsamp - burnin).", 
            call. = FALSE)
    } else if ((check <= 0) | (check > (nsamp - burnin)) | (check != round(check))) {
        stop("mcmc_control: The thinning must be a positive integer equal or less than (nsamp - burnin).", 
            call. = FALSE)
    }
    
    check <- seed
    if (!is.null(check)) {
        if ((length(check) != 1) | is.list(check)) {
            stop("mcmc_control: The seed must be an integer.", call. = FALSE)
        } else if (!is.numeric(check)) {
            stop("mcmc_control: The seed must be an integer.", call. = FALSE)
        } else if (check != round(check)) {
            stop("mcmc_control: The seed must be an integer.", call. = FALSE)
        }
    }
    
    check <- init_pars
    checknames <- c("init_alpha", "init_beta", "init_mu", "init_delta", "init_gamma", "init_nu", 
        "init_theta", "init_hidprocess")
    
    if (is.null(check)) {
        stop("mcmc_control: The init_pars has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("mcmc_control: The init_pars must be a list obtained by function mcmc_initpars.", 
            call. = FALSE)
    } else if (length(check) != 8) {
        stop("mcmc_control: The init_pars must be a list obtained by function mcmc_initpars.", 
            call. = FALSE)
    } else if (any(names(check) != checknames)) {
        stop("mcmc_control: The init_pars must be a list obtained by function mcmc_initpars.", 
            call. = FALSE)
    }
    
    check <- init_stepsize
    if (is.null(check)) {
        stop("mcmc_control: The init_stepsize has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("mcmc_control: The init_stepsize must be a positive real number.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("mcmc_control: The init_stepsize must be a positive real number.", call. = FALSE)
    } else if (check <= 0) {
        stop("mcmc_control: The init_stepsize must be a positive real number.", call. = FALSE)
    }
    
    check <- max_nsteps
    if (is.null(check)) {
        stop("mcmc_control: The max_nsteps has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("mcmc_control: The max_nsteps must be a positive integer.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("mcmc_control: The max_nsteps must be a positive integer.", call. = FALSE)
    } else if ((check <= 0) | (check != round(check))) {
        stop("mcmc_control: The max_nsteps must be a positive integer.", call. = FALSE)
    }
    
    check <- ntypes
    if (is.null(check)) {
        stop("mcmc_control: The ntypes has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("mcmc_control: The ntypes must be an integer greater than one.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("mcmc_control: The ntypes must be an integer greater than one.", call. = FALSE)
    } else if ((check <= 1) | (check != round(check))) {
        stop("mcmc_control: The ntypes must be an integer greater than one.", call. = FALSE)
    }
    
    check <- mass_diag
    if (is.null(check)) {
        stop("mcmc_control: The mass_diag has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("mcmc_control: The mass_diag must be a vector.", call. = FALSE)
    } else if (!is.null(init_pars$init_gamma) & length(check) != (3 * ntypes + 2)) {
        stop("mcmc_control: Length of mass_diag is not compatible.", call. = FALSE)
    } else if (is.null(init_pars$init_gamma) & length(check) != (3 * ntypes + 1)) {
        stop("mcmc_control: Length of mass_diag is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("mcmc_control: The entries of mass_diag must be positive reals.", call. = FALSE)
    } else if (any(check <= 0)) {
        stop("mcmc_control: The entries of mass_diag must be positive reals.", call. = FALSE)
    }
    # End of error checks for input arguments
    
    output <- list(nsamp = nsamp, burnin = burnin, thinning = thinning, seed = seed, init_pars = init_pars, 
        init_stepsize = init_stepsize, max_nsteps = max_nsteps, mass_diag = mass_diag)
    
    return(output)
}
