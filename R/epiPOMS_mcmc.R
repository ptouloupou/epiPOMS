#' MCMC-based tool for analyzing partially observed multi-strain epidemic data.
#'
#' @description
#' Runs a Bayesian data augmented Markov chain Monte Carlo (MCMC) algorithm for
#' fitting partially observed multi-strain models to longitudinal epidemic data
#' with misclassification of a population of individuals which is divided into
#' groups.
#'
#' @param
#' infoepidata A list containing the right format and general information of
#' an object of class \sQuote{epiPOMSdata}, obtained using function
#' \code{\link{info_epiPOMSdata}}.
#' @param
#' mcmccontrol A list containing the parameter values that control the MCMC
#' algorithm, obtained using function \code{\link{mcmc_control}}.
#' @param
#' priorcontrol A list containing the hyperparameters for the prior
#' distributions for the MCMC algorithm, obtained using function
#' \code{\link{prior_control}}.
#' @param
#' smallgroups A vector of positive integers, with length equal to the number
#' of groups that are smaller in terms of area (in square meters), containing
#' the indices of the small groups. These values must be between 1 and the total
#' number of groups in the study. Set \code{smallgroups = NULL} if there is no
#' difference between groups.
#' @param
#' parallel A logical value; if set to \code{TRUE}, the hidden carriage process
#' for each group is updated in parallel. If it is \code{FALSE}, then the
#' carriage processes are updated sequentially. The default is \code{TRUE}.
#' See details, below, for more information.
#' @param
#' ncores A positive integer providing the number of required cores if
#' \code{parallel} is set to \code{TRUE}. This value must be between 1 and the
#' minimum of the total number of groups in the study and the number of available
#' cores on the user's computer. By default \code{ncores} is set to the minimum of
#' the number of groups and the number of available cores.
#'
#' @details
#' Our approach to analyse partially observed longitudinal multi-strain epidemic
#' data, involves assuming that the classifications at an observation time are
#' imperfect measures of an underlying true (hidden) epidemic process. This
#' facilitates the use of partially observed epidemic models (\acronym{POMS}),
#' which provide a natural framework to analyse infection dynamics in
#' longitudinal studies where the observed data are subject to potential testing
#' error due to poor sensitivity of the diagnostic and strain procedure used.
#' For the full details see \insertCite{Touloupou2020;textual}{epiPOMS}.
#'
#' Parameter estimation in \acronym{POMS} epidemic models is done by using MCMC
#' data augmentation methods, that employs both Gibbs and Hamiltonian Monte Carlo
#' (HMC) updates. More specifically, sampling the hidden carriage states of each
#' individual is done by using a Gibbs step via the individual-Forward
#' Filtering Backward Sampling (iFFBS) algorithm by
#' \insertCite{Touloupou2019;textual}{epiPOMS}. The initial probability
#' parameters, \eqn{\nu}, and the observation parameters, \eqn{\theta_1},
#' \eqn{\theta_2}, \eqn{\theta_C}, \eqn{\theta_S} and  \eqn{\theta_P}, are
#' updated using Gibbs updates. The remaining parameters are updated jointly
#' using an HMC algorithm. For a more precise description, see
#' \code{\link{iFFBS_upd}} and \code{\link{HMC_upd}}, and
#' \code{\link{prior_control}} for the prior specification of the parameters.
#'
#' For faster run of the \code{\link{epiPOMS_mcmc}} function, the hidden
#' carriage process for each group can be updated in parallel, since we assume
#' that groups are independent of one another. This option is controlled via
#' the two arguments \code{parallel} and \code{ncores}. If  \code{parallel} is
#' set to \code{TRUE}  the number of required cores, \code{ncores}, should be
#' specified.
#'
#' @return
#' An object of class \sQuote{epiPOMSmcmc} is returned containing the following:
#' \describe{
#' \item{alphamcmc}{A \eqn{sampsize} by \eqn{n_g} matrix containing the posterior
#' samples for the \eqn{\alpha} parameters, where \eqn{sampsize} denotes the
#' size of the thinned recorded MCMC samples after the burnin period and
#' \eqn{n_g} denotes the number of different strains in the study. The
#' \eqn{i^{th}}{ith} column contains the samples for the \eqn{i^{th}}{ith}
#' \eqn{\alpha} parameter. The first column corresponds to the
#' first strain, the second column to the second strain and so on.}
#' \item{betamcmc}{A \eqn{sampsize} by \eqn{n_g} matrix containing the posterior
#' samples for the \eqn{\beta} parameters. The \eqn{i^{th}}{ith} column contains
#' the samples for the \eqn{i^{th}}{ith} \eqn{\beta} parameter. The first column
#' corresponds to the first strain, the second column to the second strain and
#' so on.}
#' \item{mumcmc}{A \eqn{sampsize} by \eqn{n_g} matrix containing the posterior
#' samples for the \eqn{\mu} parameters. The \eqn{i^{th}}{ith} column contains
#' the samples for the \eqn{i^{th}}{ith} \eqn{\mu} parameter. The first column
#' corresponds to the first strain, the second column to the second strain and
#' so on.}
#' \item{deltamcmc}{A vector of length \eqn{sampsize} containing the posterior
#' samples for parameter \eqn{\delta}.}
#' \item{gammamcmc}{A vector of length \eqn{sampsize} containing the posterior
#' samples for parameter \eqn{\gamma}, if
#' \code{smallgroups} \eqn{\neq}{≠} \code{NULL}. Otherwise a \code{NULL} value
#' is returned.}
#' \item{numcmc}{A \eqn{sampsize} by \eqn{(n_g + 1)} matrix containing the
#' posterior samples for the \eqn{\nu} parameters. The \eqn{i^{th}}{ith} column
#' contains the samples for the \eqn{i^{th}}{ith} \eqn{\nu} parameter. The first
#' column corresponds to the non-carriage state. The second and subsequent
#' columns correspond to the carriage of one of the \eqn{n_g} strains.}
#' \item{thetamcmc}{A \eqn{sampsize} by \eqn{5} matrix containing the posterior
#' samples for the \eqn{\theta} parameters. The \eqn{i^{th}}{ith} column
#' contains the samples for the \eqn{i^{th}}{ith} \eqn{\theta} parameter.}
#' \item{probstates}{A list containing \eqn{P} arrays, one for each group, with
#' the posterior probability that an individual of a group is not colonised or
#' is colonised by a specific strain over the sampling period. The dimensions of
#' each array are \eqn{(n_g + 1)} by \eqn{C^{[p]}}{C^[p]} by \eqn{T^{[p]}}{T^[p]},
#' where \eqn{C^{[p]}}{C^[p]} and \eqn{T^{[p]}}{T^[p]} denote the total number of
#' individuals and the time of the last observation in group
#' \eqn{p = 1, 2, \ldots, P}, respectively.}
#' \item{sumstates}{A \eqn{sampsize} by \eqn{(n_g + 1)} matrix containing the
#' total number of individuals in the population having augmented state of
#' carriage equal to \eqn{s} over the entire study period, where
#' \eqn{s = 0, 1, \ldots, n_g}.}
#' \item{statetrans}{A \eqn{sampsize} by \eqn{5} matrix containing the total
#' number of state-specific transitions in the augmented carriage process of all
#' individuals. The \eqn{i^{th}}{ith} column contains the total number for the
#' \eqn{i^{th}}{ith} state-specific transition. The transitions are (in order):
#' from \eqn{r} to \eqn{s} where \eqn{r \neq s \neq 0}{r ≠ s ≠ 0}, from \eqn{0}
#' to \eqn{s} where \eqn{s \neq 0}{s ≠ 0}, from \eqn{r} to \eqn{0} where
#' \eqn{r \neq 0}{r ≠ 0}, from \eqn{0} to \eqn{0} and from \eqn{r} to \eqn{r}
#' where \eqn{r \neq 0}{r ≠ 0}.}
#' \item{stepsize}{A real value specifying the leapfrog stepsize used in HMC
#' method after the burnin period.}
#' }
#'
#' @author
#' Panayiota Touloupou
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{prior_control}} for specifying the hyperparameters of the prior
#' distributions, \code{\link{mcmc_control}} for setting the parameter values that
#' control the MCMC algorithm, \code{\link{info_epiPOMSdata}} for extracting
#' information of an object of class \sQuote{epiPOMSdata},
#' \code{\link{plot.epiPOMSmcmc}} and \code{\link{summary.epiPOMSmcmc}} for
#' plotting and displaying summary information of an \sQuote{epiPOMSmcmc}
#' object.
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
#' \dontrun{
#' # Note: This will take a few minutes to run
#' # E. coli O157:H7 data analysis
#' # Load E. coli O157:H7 data
#' set.seed(50)
#' data(Ecoliepidata)
#' smallgroups <- c(6:10, 14:20)
#'
#' # Format the data
#' epidata <- as_epiPOMSdata(Ecoliepidata)
#'
#' # Modify and extract information
#' infoepidata <- info_epiPOMSdata(epidata)
#'
#' # Set hyperparameters for the prior distributions
#' priorcontrol <- prior_control(ntypes = infoepidata$ntypes)
#'
#' # Generate initial values for MCMC algorithm
#' mcmcinitpars <- mcmc_initpars(priorcontrol = priorcontrol, epidata =
#'      epidata, smallgroups = smallgroups)
#'
#' # Set control parameters for MCMC algorithm
#' mcmccontrol <- mcmc_control(nsamp = 250, burnin = 50, init_pars =
#'      mcmcinitpars, ntypes = infoepidata$ntypes)
#'
#' # Run MCMC algorithm on the E. coli O157:H7 data
#' # Note: Not enough iterations for any real inference
#' tic = Sys.time()
#' examplemcmc <- epiPOMS_mcmc(infoepidata = infoepidata, mcmccontrol =
#'      mcmccontrol, priorcontrol = priorcontrol, smallgroups = smallgroups,
#'      parallel = TRUE)
#' toc = Sys.time() - tic
#' toc
#'}
#'
#' \dontrun{
#' # Note: This will take a few minutes to run
#' # E. coli O157:H7 data analysis using a simpler model where parameter delta
#' # is not included in the model
#' # Load E. coli O157:H7 data
#' set.seed(50)
#' data(Ecoliepidata)
#' smallgroups <- NULL # we assume no difference between small and big pens
#'
#' # Format the data
#' epidata <- as_epiPOMSdata(Ecoliepidata)
#'
#' # Modify and extract information
#' infoepidata <- info_epiPOMSdata(epidata)
#'
#' # Set hyperparameters for the prior distributions
#' priorcontrol <- prior_control(ntypes = infoepidata$ntypes, gammaprior = NULL)
#'
#' # Generate initial values for MCMC algorithm
#' mcmcinitpars <- mcmc_initpars(priorcontrol = priorcontrol, epidata =
#'      epidata, smallgroups = smallgroups)
#'
#' # Set control parameters for MCMC algorithm
#' mcmccontrol <- mcmc_control(nsamp = 250, burnin = 50, init_pars =
#'      mcmcinitpars, ntypes = infoepidata$ntypes)
#'
#' # Run MCMC algorithm on the E. coli O157:H7 data
#' # Note: Not enough iterations for any real inference
#' tic = Sys.time()
#' examplemcmc <- epiPOMS_mcmc(infoepidata = infoepidata, mcmccontrol =
#'      mcmccontrol, priorcontrol = priorcontrol, smallgroups = smallgroups,
#'      parallel = TRUE)
#' toc = Sys.time() - tic
#' toc
#' }
#'
#' @export
epiPOMS_mcmc <- function(infoepidata, mcmccontrol, priorcontrol, smallgroups, parallel = TRUE, 
    ncores = min(infoepidata$ngroups, getOption("cl.cores", detectCores()))) {
    
    # Error checks for input arguments
    check <- infoepidata
    checknames <- c("ntypes", "obserdata", "ngroups", "ninds", "tmax", "indNA", "tmaxNA")
    if (is.null(check)) {
        stop("epiPOMS_mcmc: The infoepidata has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("epiPOMS_mcmc: The infoepidata must be a list obtained by function info_epiPOMSdata.", 
            call. = FALSE)
    } else if (length(check) != 7) {
        stop("epiPOMS_mcmc: The infoepidata must be a list obtained by function info_epiPOMSdata.", 
            call. = FALSE)
    } else if (any(names(check) != checknames)) {
        stop("epiPOMS_mcmc: The infoepidata must be a list obtained by function info_epiPOMSdata.", 
            call. = FALSE)
    }
    
    check <- mcmccontrol
    checknames <- c("nsamp", "burnin", "thinning", "seed", "init_pars", "init_stepsize", "max_nsteps", 
        "mass_diag")
    if (is.null(check)) {
        stop("epiPOMS_mcmc: The mcmccontrol has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("epiPOMS_mcmc: The mcmccontrol must be a list obtained by function mcmc_control.", 
            call. = FALSE)
    } else if (length(check) != 8) {
        stop("epiPOMS_mcmc: The mcmccontrol must be a list obtained by function mcmc_control.", 
            call. = FALSE)
    } else if (any(names(check) != checknames)) {
        stop("epiPOMS_mcmc: The mcmccontrol must be a list obtained by function mcmc_control.", 
            call. = FALSE)
    }
    
    check <- priorcontrol
    checknames <- c("alphaprior", "betaprior", "muprior", "deltaprior", "gammaprior", "nuprior", 
        "theta12prior", "thetaCSprior", "thetaPprior")
    if (is.null(check)) {
        stop("epiPOMS_mcmc: The priorcontrol has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("epiPOMS_mcmc: The priorcontrol must be a list obtained by function prior_control.", 
            call. = FALSE)
    } else if (length(check) != 9) {
        stop("epiPOMS_mcmc: The priorcontrol must be a list obtained by function prior_control.", 
            call. = FALSE)
    } else if (any(names(check) != checknames)) {
        stop("epiPOMS_mcmc: The priorcontrol must be a list obtained by function prior_control.", 
            call. = FALSE)
    }
    
    if ((is.null(mcmccontrol$init_pars$init_gamma)) != (is.null(priorcontrol$gammaprior))) {
        stop("epiPOMS_mcmc: The init_gamma at mcmccontrol$init_pars must be NULL if gammaprior at priorcontrol is NULL, and via versa.", 
            call. = FALSE)
    }
    
    check <- smallgroups
    if (!is.null(check)) {
        if ((!is.vector(check)) | is.list(check)) {
            stop("epiPOMS_mcmc: The smallgroups must be a vector.", call. = FALSE)
        } else if (length(check) > infoepidata$ngroups) {
            stop("epiPOMS_mcmc: Length of smallgroups is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("epiPOMS_mcmc: The entries of smallgroups must be integers between 1 and the total number of groups.", 
                call. = FALSE)
        } else if (any(check < 1) | any(check > infoepidata$ngroups) | any(check != round(check))) {
            stop("epiPOMS_mcmc: The entries of smallgroups must be integers between 1 and the total number of groups.", 
                call. = FALSE)
        }
    }
    
    if ((is.null(smallgroups)) != (is.null(priorcontrol$gammaprior))) {
        stop("epiPOMS_mcmc: The smallgroups must be NULL if gammaprior at priorcontrol is NULL, and via versa.", 
            call. = FALSE)
    }
    
    check <- parallel
    if (is.null(check)) {
        stop("epiPOMS_mcmc: The parallel has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("epiPOMS_mcmc: The parallel must be logical, either TRUE or FALSE.", call. = FALSE)
    } else if (!is.logical(check) | is.na(check)) {
        stop("epiPOMS_mcmc: The parallel must be logical, either TRUE or FALSE.", call. = FALSE)
    }
    
    if (parallel == TRUE) {
        check <- ncores
        if (is.null(check)) {
            stop("epiPOMS_mcmc: The ncores has to be specified.", call. = FALSE)
        } else if ((length(check) != 1) | is.list(check)) {
            stop("epiPOMS_mcmc: The ncores must be an integer between 1 and the minimum of the number of groups and the number of available cores.", 
                call. = FALSE)
        } else if (!is.numeric(check)) {
            stop("epiPOMS_mcmc: The ncores must be an integer between 1 and the minimum of the number of groups and the number of available cores.", 
                call. = FALSE)
        } else if ((check < 1) | (check != round(check)) | (check > min(infoepidata$ngroups, getOption("cl.cores", 
            detectCores())))) {
            stop("epiPOMS_mcmc: The ncores must be an integer between 1 and the minimum of the number of groups and the number of available cores.", 
                call. = FALSE)
        }
        # ncores <- ncores - 1
        ncores <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(ncores)
        on.exit(parallel::stopCluster(ncores), add = TRUE)
        
    }
    # End of error checks for input arguments
    
    if (!is.null(mcmccontrol$seed)) {
        set.seed(mcmccontrol$seed)
    }
    
    sampsize <- rnd((mcmccontrol$nsamp - mcmccontrol$burnin)/(mcmccontrol$thinning))
    
    alphamcmc <- matrix(NA, nrow = sampsize, ncol = infoepidata$ntypes)
    betamcmc <- matrix(NA, nrow = sampsize, ncol = infoepidata$ntypes)
    mumcmc <- matrix(NA, nrow = sampsize, ncol = infoepidata$ntypes)
    deltamcmc <- rep(NA, sampsize)
    if (!is.null(smallgroups)) {
        gammamcmc <- rep(NA, sampsize)
    } else {
        gammamcmc <- NULL
    }
    numcmc <- matrix(NA, nrow = sampsize, ncol = (infoepidata$ntypes + 1))
    thetamcmc <- matrix(NA, nrow = sampsize, ncol = 5)
    
    
    logcuralpha <- log(mcmccontrol$init_pars$init_alpha)
    logcurbeta <- log(mcmccontrol$init_pars$init_beta)
    logcurmu <- log(mcmccontrol$init_pars$init_mu)
    logcurdelta <- log(mcmccontrol$init_pars$init_delta)
    if (!is.null(smallgroups)) {
        logcurgamma <- log(mcmccontrol$init_pars$init_gamma)
    } else {
        logcurgamma <- mcmccontrol$init_pars$init_gamma
    }
    curnu <- mcmccontrol$init_pars$init_nu
    curobserpars <- mcmccontrol$init_pars$init_theta
    curhidproc <- mcmccontrol$init_pars$init_hidprocess
    curclassmat <- classprob(infoepidata$ntypes, curobserpars)
    
    
    probstates <- list()
    for (kk in 1:infoepidata$ngroups) {
        probstates[[kk]] <- array(0, dim = c((infoepidata$ntypes + 1), infoepidata$ninds[kk], infoepidata$tmax[kk]))
        
    }
    
    mstatetrans <- matrix(NA, nrow = sampsize, ncol = 5)
    msumstates <- matrix(NA, nrow = sampsize, ncol = (infoepidata$ntypes + 1))
    
    posNA <- which(is.na(unlist(infoepidata$obserdata$test1)))
    ithinning <- 0
    epsilon <- mcmccontrol$init_stepsize
    accforepsilon <- rep(NA, mcmccontrol$nsamp - mcmccontrol$burnin)
    
    for (i in 1:mcmccontrol$nsamp) {
        ### Sample hidden states
        
        if (parallel == TRUE) {
            out <- list()
            out <- foreach(kk = 1:infoepidata$ngroups) %dopar% {
                gamma_new <- ifelse(kk %in% smallgroups, exp(logcurgamma), 1)
                curhidprocgr <- curhidproc[[kk]]
                for (jj in 1:infoepidata$ninds[kk]) {
                  curhidprocgr[jj, ] <- iFFBS_upd(ntypes = infoepidata$ntypes, nindsgr = infoepidata$ninds[kk], 
                    ind = jj, tmaxgr = infoepidata$tmax[kk], tmaxind = infoepidata$tmaxNA[[kk]][jj], 
                    curhidprocgr = curhidprocgr, obsertest1ind = infoepidata$obserdata$test1[[kk]][jj, 
                      ], obsertest2ind = infoepidata$obserdata$test2[[kk]][jj, ], curalpha = exp(logcuralpha), 
                    curbeta = exp(logcurbeta), curmu = exp(logcurmu), curdelta = exp(logcurdelta), 
                    curgamma = gamma_new, curnu = curnu, curclassmat = curclassmat)
                }
                return(curhidprocgr)
            }
            curhidproc <- out
        } else {
            for (kk in 1:infoepidata$ngroups) {
                gamma_new <- ifelse(kk %in% smallgroups, exp(logcurgamma), 1)
                for (jj in 1:infoepidata$ninds[kk]) {
                  curhidproc[[kk]][jj, ] <- iFFBS_upd(ntypes = infoepidata$ntypes, nindsgr = infoepidata$ninds[kk], 
                    ind = jj, tmaxgr = infoepidata$tmax[kk], tmaxind = infoepidata$tmaxNA[[kk]][jj], 
                    curhidprocgr = curhidproc[[kk]], obsertest1ind = infoepidata$obserdata$test1[[kk]][jj, 
                      ], obsertest2ind = infoepidata$obserdata$test2[[kk]][jj, ], curalpha = exp(logcuralpha), 
                    curbeta = exp(logcurbeta), curmu = exp(logcurmu), curdelta = exp(logcurdelta), 
                    curgamma = gamma_new, curnu = curnu, curclassmat = curclassmat)
                }
            }
        }
        
        ### Hamiltonian Monte Carlo updates
        HMCoutput <- HMC_upd(ntypes = infoepidata$ntypes, ngroups = infoepidata$ngroups, ninds = infoepidata$ninds, 
            tmax = infoepidata$tmax, logcuralpha = logcuralpha, logcurbeta = logcurbeta, logcurmu = logcurmu, 
            logcurdelta = logcurdelta, logcurgamma = logcurgamma, curhidproc = curhidproc, stepsize = epsilon, 
            nsteps = mcmccontrol$max_nsteps, mass_diag = mcmccontrol$mass_diag, priorcontrol = priorcontrol, 
            smallgroups = smallgroups)
        logcuralpha <- HMCoutput$logcuralpha
        logcurbeta <- HMCoutput$logcurbeta
        logcurmu <- HMCoutput$logcurmu
        logcurdelta <- HMCoutput$logcurdelta
        logcurgamma <- HMCoutput$logcurgamma
        
        ### Gibbs step for updating the initial probability parameters
        Ininf <- rep(0, (infoepidata$ntypes + 1))
        zInfec <- as.vector(unlist(lapply(curhidproc, function(x) x[, 1])))
        tI <- table(zInfec)
        pos <- as.numeric(names(tI)) + 1
        Ininf[pos] <- tI
        curnu <- MCMCpack::rdirichlet(n = 1, alpha = Ininf + priorcontrol$nuprior)
        
        ### Gibbs step for updating the observation parameters
        obserdatatest1help <- unlist(infoepidata$obserdata$test1)[-posNA]
        obserdatatest2help <- unlist(infoepidata$obserdata$test2)[-posNA]
        hidprochelp <- unlist(curhidproc)[-posNA]
        obsparshelp <- .Call("updateobspars", as.integer(length(obserdatatest1help)), as.integer(obserdatatest1help), 
            as.integer(obserdatatest2help), as.integer(hidprochelp), as.integer(infoepidata$ntypes))
        curobserpars[1] <- rbeta(n = 1, shape1 = priorcontrol$theta12prior[1] + obsparshelp[1], 
            shape2 = priorcontrol$theta12prior[2] + obsparshelp[2])
        curobserpars[2] <- rbeta(n = 1, shape1 = priorcontrol$theta12prior[1] + obsparshelp[3], 
            shape2 = priorcontrol$theta12prior[2] + obsparshelp[4])
        curobserpars[3:4] <- MCMCpack::rdirichlet(n = 1, alpha = obsparshelp[5:7] + priorcontrol$thetaCSprior)[1:2]
        curobserpars[5] <- rbeta(n = 1, shape1 = priorcontrol$thetaPprior[1] + obsparshelp[8], 
            shape2 = priorcontrol$thetaPprior[2] + obsparshelp[9])
        
        curclassmat <- classprob(infoepidata$ntypes, curobserpars)
        
        
        # adjust stepsize of HMC in burnin period
        accforepsilon[i] <- logcurdelta
        if (((i%%100) == 0) & (i <= mcmccontrol$burnin)) {
            acc <- length(which(diff(accforepsilon[(i - 99):i]) != 0))/99
            y <- 1 + 1000 * (acc - 0.7) * (acc - 0.7) * (acc - 0.7)
            if (y < 0.9) {
                epsilon <- 0.9 * epsilon
            } else if (y > 1.1) {
                epsilon <- 1.1 * epsilon
            }
        }
        
        # recorded MCMC samples
        if ((i > mcmccontrol$burnin) & ((i%%mcmccontrol$thinning) == 0)) {
            ithinning <- ithinning + 1
            
            alphamcmc[ithinning, ] <- exp(logcuralpha)
            betamcmc[ithinning, ] <- exp(logcurbeta)
            mumcmc[ithinning, ] <- exp(logcurmu)
            deltamcmc[ithinning] <- exp(logcurdelta)
            if (!is.null(smallgroups)) {
                gammamcmc[ithinning] <- exp(logcurgamma)
                
            }
            numcmc[ithinning, ] <- curnu
            thetamcmc[ithinning, ] <- curobserpars
            
            statetranshelp <- rep(0, 5)
            for (kk in 1:infoepidata$ngroups) {
                hidproc_vec <- as.vector(curhidproc[[kk]])
                hidproc_vec[which(is.na(hidproc_vec))] <- -1
                statetranshelp <- statetranshelp + .Call("statetrans", as.integer(infoepidata$tmax[kk]), 
                  as.integer(infoepidata$ninds[kk]), as.integer(hidproc_vec))
            }
            mstatetrans[ithinning, ] <- statetranshelp
            
            hidproc_vec <- as.vector(unlist(curhidproc))
            hidproc_vec <- hidproc_vec[-which(is.na(hidproc_vec))]
            msumstates[ithinning, ] <- .Call("sumstates", as.integer(hidproc_vec), as.integer(length(hidproc_vec)), 
                as.integer(infoepidata$ntypes + 1))
            
            for (kk in 1:infoepidata$ngroups) {
                for (ii in 1:infoepidata$tmax[kk]) {
                  for (jj in 1:infoepidata$ninds[kk]) {
                    for (ts in 0:infoepidata$ntypes) {
                      if (!is.na(curhidproc[[kk]][jj, ii]) & curhidproc[[kk]][jj, ii] == ts) {
                        probstates[[kk]][ts + 1, jj, ii] <- probstates[[kk]][ts + 1, jj, ii] + 
                          1
                      }
                    }
                  }
                }
                
            }
            
        }
        
    }
    
    
    output <- list(alphamcmc = alphamcmc, betamcmc = betamcmc, mumcmc = mumcmc, deltamcmc = deltamcmc, 
        gammamcmc = gammamcmc, numcmc = numcmc, thetamcmc = thetamcmc, probstates = lapply(probstates, 
            function(x) x/sampsize), sumstates = msumstates, statetrans = mstatetrans, stepsize = epsilon)
    class(output) <- "epiPOMSmcmc"
    return(output)
}


rnd <- function(x) trunc(x + sign(x) * 0.5)





