#' Update the transmission parameters using an HMC algorithm.
#'
#' @description
#' This function updates the transmission parameters using an Hamiltonian Monte
#' Carlo (HMC) algorithm. This function is only used in conjunction with
#' \code{\link{epiPOMS_mcmc}}.
#'
#' @inheritParams hidstate_sim
#' @param
#' logcuralpha A vector of real numbers, one for each strain, containing the
#' current values of the logarithm of \eqn{\alpha} parameters, i.e. the
#' strain-specific external colonisation rates. This is a vector of \eqn{n_g}
#' values. The first element corresponds to the first strain, the second element
#' to the second strain and so on.
#' @param
#' logcurbeta A vector of real numbers, one for each strain, containing the
#' current values of the logarithm of \eqn{\beta} parameters, i.e. the
#' strain-specific within-group colonisation rates. This is a vector of \eqn{n_g}
#' values. The first element corresponds to the first strain, the second element
#' to the second strain and so on.
#' @param
#' logcurmu A vector of real numbers, one for each strain, containing the
#' current values of the logarithm of \eqn{\mu} parameters, i.e. the
#' strain-specific clearance rates. This is a vector of \eqn{n_g} values.
#' The first element corresponds to the first strain, the second element to
#' the second strain and so on.
#' @param
#' logcurdelta A real number corresponding to the current value of the logarithm
#' of \eqn{\delta} parameter, i.e. the relative colonisation rate in a
#' carrier versus non-carrier individual.
#' @param
#' logcurgamma A real number corresponding to the current value of the logarithm
#' of \eqn{\gamma} parameter, i.e. the relative colonisation rate in
#' smaller versus bigger groups in terms of area (in square meters). If
#' \code{smallgroups = NULL}, i.e there is no difference between groups, set
#' \code{logcurgamma = NULL}.
#' @param
#' curhidproc A list containing \eqn{P} matrices (one for each group), with
#' \code{ninds}\eqn{[p]} rows and \code{tmax}\eqn{[p]} columns
#' for each \eqn{p = 1, 2, \ldots, P}. Each matrix contains the current carriage
#' states of all individuals within the group at each time point over their
#' group observation period. All elements must be non-negative integers, ranging
#' from 0 to \eqn{n_g}, or \code{NA}. \code{NA} entries corresponds to missing
#' values, for example, individual dropouts.
#' @param
#' stepsize A positive real value specifying the leapfrog stepsize in HMC
#' method.
#' @param
#' nsteps A positive integer corresponding to the number of leapfrog steps in
#' HMC method.
#' @param
#' mass_diag A vector of positive numbers providing the diagonal of the mass
#' matrix, with length equal to the number of parameters to be updated jointly
#' using the HMC algorithm, i.e. \eqn{(3 \times n_g + 2)}{(3xn_g + 2)} if
#' \code{smallgroups} \eqn{\neq}{≠} \code{NULL} or
#' \eqn{(3 \times n_g + 1)}{(3xn_g + 1)} if \code{smallgroups = NULL}. See
#' details, below, for more information.
#' @param
#' priorcontrol A list containing the hyperparameters for the prior
#' distributions for the MCMC algorithm, obtained using function
#' \code{\link{prior_control}}.
#'
#' @details
#' Sampling from the posterior distribution of the parameters of partially
#' observed multi-strain epidemic models is performed within a Bayesian
#' framework using an MCMC algorithm (see
#' \insertCite{Touloupou2020;textual}{epiPOMS} for more details). Some of the
#' full conditionals are not given in closed form and therefore for these we
#' resort to HMC \insertCite{Neal2011}{epiPOMS}. More specifically, the following
#' parameters are updated jointly, in order, with HMC: the strain-specific
#' external colonisation rates (\eqn{\alpha}), the within-group colonisation
#' rates (\eqn{\beta}), the clearance rates (\eqn{\mu}), the relative
#' colonisation rate in a carrier versus non-carrier individual (\eqn{\delta})
#' and, if \code{smallgroups} \eqn{\neq}{≠} \code{NULL}, the relative colonisation
#' rate in smaller versus bigger groups in terms of area (\eqn{\gamma}). An
#' important constraint on these model parameters is that they must be positive
#' quantities. Additionally, some of the rates may become very small and therefore
#' we follow the common practice of using a logarithmic transformation for these
#' parameters.
#'
#' In the sampling process of HMC, leapfrog integration is used to
#' approximately solve Hamilton's equations, and the integration is required
#' to set the number of discrete time steps (\code{nsteps}) and the integration
#' stepsize (\code{stepsize}). Moreover, the specification of a mass matrix is
#' required which should be symmetric and positive-definite, and is often a scalar
#' multiple of the identity matrix. Here, we assume that the mass matrix is
#' diagonal, with diagonal elements being the elements of \code{mass_diag}.
#'
#' @return
#' A list is returned with the HMC output of the following updated model
#' parameters:
#' \describe{
#' \item{logcuralpha}{A vector with the updated values of the logarithm of
#' \eqn{\alpha} parameters.}
#' \item{logcurbeta}{A vector with the updated values of the logarithm of
#' \eqn{\beta} parameters.}
#' \item{logcurmu}{A vector with the updated values of the logarithm of
#' \eqn{\mu} parameters.}
#' \item{logcurdelta}{The updated value of the logarithm of
#' \eqn{\delta} parameter.}
#' \item{logcurgamma}{The updated value of the logarithm of
#' \eqn{\gamma} parameter, if \code{smallgroups} \eqn{\neq}{≠} \code{NULL}.
#' Otherwise a \code{NULL} value is returned.}
#' }
#'
#' @author
#' Panayiota Touloupou
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for performing inference in partially observed
#' multi-strain epidemic  models.
#'
#' @export
HMC_upd <- function(ntypes, ngroups, ninds, tmax, logcuralpha, logcurbeta, logcurmu, logcurdelta, 
    logcurgamma, curhidproc, stepsize, nsteps, mass_diag, priorcontrol, smallgroups) {
    
    if (!(is.null(smallgroups))) {
        q <- c(logcuralpha, logcurbeta, logcurmu, logcurdelta, logcurgamma)
        p <- rnorm(length(q), 0, 1)
        p <- p * sqrt(mass_diag)
        curp <- p
        
        alphaprior_rates <- rep(x = priorcontrol$alphaprior, times = ntypes)
        betaprior_rates <- rep(x = priorcontrol$betaprior, times = ntypes)
        muprior_rates <- rep(x = priorcontrol$muprior, times = ntypes)
        
        logpost1 <- 0
        for (kk in 1:ngroups) {
            indicatorsmall <- ifelse(kk %in% smallgroups, 1, 0)
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            logpost1 <- logpost1 + .Call("logpost", as.double(exp(c(0, logcuralpha))), as.double(exp(c(0, 
                logcurbeta))), as.double(exp(c(0, logcurmu))), as.double(exp(logcurdelta)), as.double(exp(logcurgamma)), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)))
        }
        
        
        logpost1 <- logpost1 + sum(dexp(exp(logcuralpha), rate = alphaprior_rates, log = TRUE) + 
            log(exp(logcuralpha)))
        logpost1 <- logpost1 + sum(dexp(exp(logcurbeta), rate = betaprior_rates, log = TRUE) + 
            log(exp(logcurbeta)))
        logpost1 <- logpost1 + sum(dexp(exp(logcurmu), rate = muprior_rates, log = TRUE) + log(exp(logcurmu)))
        logpost1 <- logpost1 + dexp(exp(logcurdelta), rate = priorcontrol$deltaprior, log = TRUE) + 
            log(exp(logcurdelta))
        logpost1 <- logpost1 + dexp(exp(logcurgamma), rate = priorcontrol$gammaprior, log = TRUE) + 
            log(exp(logcurgamma))
        
        CurrentH <- -logpost1 + 0.5 * sum((curp^2)/(mass_diag))
        
        grad1 <- rep(0, 3 * ntypes + 2)
        for (kk in 1:ngroups) {
            indicatorsmall <- ifelse(kk %in% smallgroups, 1, 0)
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            last <- ifelse(kk == ngroups, 1, 0)
            grad1 <- grad1 + .Call("gradient", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * ntypes)])), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)), as.double(c(priorcontrol$alphaprior, priorcontrol$betaprior, 
                  priorcontrol$muprior, priorcontrol$deltaprior, priorcontrol$gammaprior)), as.integer(last))
        }
        
        p <- p + stepsize * grad1/2
        nsteps <- ceiling(runif(1) * nsteps) - 1
        
        for (i in 1:nsteps) {
            q <- q + stepsize * p/(mass_diag)
            
            if (i != nsteps) {
                
                grad1 <- rep(0, 3 * ntypes + 2)
                for (kk in 1:ngroups) {
                  indicatorsmall <- ifelse(kk %in% smallgroups, 1, 0)
                  curhidproc_vec <- as.vector(curhidproc[[kk]])
                  last <- ifelse(kk == ngroups, 1, 0)
                  grad1 <- grad1 + .Call("gradient", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                    q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                    ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * 
                    ntypes)])), as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), 
                    as.integer(curhidproc_vec), as.integer(ntypes + 1), as.integer(rep(0, ntypes + 
                      1)), as.integer(matrix(0, ntypes + 1, ntypes + 1)), as.double(c(priorcontrol$alphaprior, 
                      priorcontrol$betaprior, priorcontrol$muprior, priorcontrol$deltaprior, priorcontrol$gammaprior)), 
                    as.integer(last))
                }
                
                p <- p + stepsize * grad1
                
            }
        }
        
        grad1 <- rep(0, 3 * ntypes + 2)
        for (kk in 1:ngroups) {
            indicatorsmall <- ifelse(kk %in% smallgroups, 1, 0)
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            last <- ifelse(kk == ngroups, 1, 0)
            grad1 <- grad1 + .Call("gradient", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * ntypes)])), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)), as.double(c(priorcontrol$alphaprior, priorcontrol$betaprior, 
                  priorcontrol$muprior, priorcontrol$deltaprior, priorcontrol$gammaprior)), as.integer(last))
        }
        p <- p + stepsize * grad1/2
        p <- -p
        
        logpost2 <- 0
        for (kk in 1:ngroups) {
            indicatorsmall <- ifelse(kk %in% smallgroups, 1, 0)
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            logpost2 <- logpost2 + .Call("logpost", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * ntypes)])), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)))
        }
        
        
        logpost2 <- logpost2 + sum(dexp(exp(q[1:(1 * ntypes)]), rate = alphaprior_rates, log = TRUE) + 
            log(exp(q[1:(1 * ntypes)])))
        logpost2 <- logpost2 + sum(dexp(exp(q[(1 * ntypes + 1):(2 * ntypes)]), rate = betaprior_rates, 
            log = TRUE) + log(exp(q[(1 * ntypes + 1):(2 * ntypes)])))
        logpost2 <- logpost2 + sum(dexp(exp(q[(2 * ntypes + 1):(3 * ntypes)]), rate = muprior_rates, 
            log = TRUE) + log(exp(q[(2 * ntypes + 1):(3 * ntypes)])))
        logpost2 <- logpost2 + dexp(exp(q[(1 + 3 * ntypes)]), rate = priorcontrol$deltaprior, log = TRUE) + 
            log(exp(q[(1 + 3 * ntypes)]))
        logpost2 <- logpost2 + dexp(exp(q[(2 + 3 * ntypes)]), rate = priorcontrol$gammaprior, log = TRUE) + 
            log(exp(q[(2 + 3 * ntypes)]))
        
        ProposedH <- -logpost2 + 0.5 * sum((p^2)/(mass_diag))
        
        ap <- -ProposedH + CurrentH
        alpha <- min(1, exp(ap))
        u <- runif(1)
        if (is.na(alpha)) {
            alpha <- 0
        }
        if (u < alpha) {
            output <- list(logcuralpha = q[1:(1 * ntypes)], logcurbeta = q[(1 * ntypes + 1):(2 * 
                ntypes)], logcurmu = q[(2 * ntypes + 1):(3 * ntypes)], logcurdelta = q[(1 + 3 * 
                ntypes)], logcurgamma = q[(2 + 3 * ntypes)])
            return(output)
        } else {
            output <- list(logcuralpha = logcuralpha, logcurbeta = logcurbeta, logcurmu = logcurmu, 
                logcurdelta = logcurdelta, logcurgamma = logcurgamma)
            return(output)
        }
    } else {
        logcurgamma <- 0
        priorcontrol$gammaprior <- 0
        qprop <- c(logcuralpha, logcurbeta, logcurmu, logcurdelta)
        p <- rnorm(length(qprop), 0, 1)
        p <- p * sqrt(mass_diag)
        curp <- p
        
        alphaprior_rates <- rep(x = priorcontrol$alphaprior, times = ntypes)
        betaprior_rates <- rep(x = priorcontrol$betaprior, times = ntypes)
        muprior_rates <- rep(x = priorcontrol$muprior, times = ntypes)
        
        logpost1 <- 0
        for (kk in 1:ngroups) {
            indicatorsmall <- 0
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            logpost1 <- logpost1 + .Call("logpost", as.double(exp(c(0, logcuralpha))), as.double(exp(c(0, 
                logcurbeta))), as.double(exp(c(0, logcurmu))), as.double(exp(logcurdelta)), as.double(exp(logcurgamma)), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)))
        }
        
        
        logpost1 <- logpost1 + sum(dexp(exp(logcuralpha), rate = alphaprior_rates, log = TRUE) + 
            log(exp(logcuralpha)))
        logpost1 <- logpost1 + sum(dexp(exp(logcurbeta), rate = betaprior_rates, log = TRUE) + 
            log(exp(logcurbeta)))
        logpost1 <- logpost1 + sum(dexp(exp(logcurmu), rate = muprior_rates, log = TRUE) + log(exp(logcurmu)))
        logpost1 <- logpost1 + dexp(exp(logcurdelta), rate = priorcontrol$deltaprior, log = TRUE) + 
            log(exp(logcurdelta))
        
        CurrentH <- -logpost1 + 0.5 * sum((curp^2)/(mass_diag))
        
        q <- c(qprop[1:(3 * ntypes + 1)], logcurgamma)
        grad1 <- rep(0, 3 * ntypes + 2)
        for (kk in 1:ngroups) {
            indicatorsmall <- 0
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            last <- ifelse(kk == ngroups, 1, 0)
            grad1 <- grad1 + .Call("gradient", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * ntypes)])), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)), as.double(c(priorcontrol$alphaprior, priorcontrol$betaprior, 
                  priorcontrol$muprior, priorcontrol$deltaprior, priorcontrol$gammaprior)), as.integer(last))
        }
        grad1 <- grad1[-(2 + 3 * ntypes)]
        
        p <- p + stepsize * grad1/2
        nsteps <- ceiling(runif(1) * nsteps) - 1
        
        for (i in 1:nsteps) {
            qprop <- qprop + stepsize * p/(mass_diag)
            
            if (i != nsteps) {
                
                q <- c(qprop[1:(3 * ntypes + 1)], logcurgamma)
                grad1 <- rep(0, 3 * ntypes + 2)
                for (kk in 1:ngroups) {
                  indicatorsmall <- 0
                  curhidproc_vec <- as.vector(curhidproc[[kk]])
                  last <- ifelse(kk == ngroups, 1, 0)
                  grad1 <- grad1 + .Call("gradient", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                    q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                    ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * 
                    ntypes)])), as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), 
                    as.integer(curhidproc_vec), as.integer(ntypes + 1), as.integer(rep(0, ntypes + 
                      1)), as.integer(matrix(0, ntypes + 1, ntypes + 1)), as.double(c(priorcontrol$alphaprior, 
                      priorcontrol$betaprior, priorcontrol$muprior, priorcontrol$deltaprior, priorcontrol$gammaprior)), 
                    as.integer(last))
                }
                grad1 <- grad1[-(2 + 3 * ntypes)]
                
                p <- p + stepsize * grad1
                
            }
        }
        
        q <- c(qprop[1:(3 * ntypes + 1)], logcurgamma)
        grad1 <- rep(0, 3 * ntypes + 2)
        for (kk in 1:ngroups) {
            indicatorsmall <- 0
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            last <- ifelse(kk == ngroups, 1, 0)
            grad1 <- grad1 + .Call("gradient", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * ntypes)])), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)), as.double(c(priorcontrol$alphaprior, priorcontrol$betaprior, 
                  priorcontrol$muprior, priorcontrol$deltaprior, priorcontrol$gammaprior)), as.integer(last))
        }
        grad1 <- grad1[-(2 + 3 * ntypes)]
        
        p <- p + stepsize * grad1/2
        p <- -p
        
        logpost2 <- 0
        for (kk in 1:ngroups) {
            indicatorsmall <- 0
            curhidproc_vec <- as.vector(curhidproc[[kk]])
            logpost2 <- logpost2 + .Call("logpost", as.double(exp(c(0, q[1:(1 * ntypes)]))), as.double(exp(c(0, 
                q[(1 * ntypes + 1):(2 * ntypes)]))), as.double(exp(c(0, q[(2 * ntypes + 1):(3 * 
                ntypes)]))), as.double(exp(q[(1 + 3 * ntypes)])), as.double(exp(q[(2 + 3 * ntypes)])), 
                as.integer(ninds[kk]), as.integer(tmax[kk]), as.integer(indicatorsmall), as.integer(curhidproc_vec), 
                as.integer(ntypes + 1), as.integer(rep(0, ntypes + 1)), as.integer(matrix(0, ntypes + 
                  1, ntypes + 1)))
        }
        
        
        logpost2 <- logpost2 + sum(dexp(exp(q[1:(1 * ntypes)]), rate = alphaprior_rates, log = TRUE) + 
            log(exp(q[1:(1 * ntypes)])))
        logpost2 <- logpost2 + sum(dexp(exp(q[(1 * ntypes + 1):(2 * ntypes)]), rate = betaprior_rates, 
            log = TRUE) + log(exp(q[(1 * ntypes + 1):(2 * ntypes)])))
        logpost2 <- logpost2 + sum(dexp(exp(q[(2 * ntypes + 1):(3 * ntypes)]), rate = muprior_rates, 
            log = TRUE) + log(exp(q[(2 * ntypes + 1):(3 * ntypes)])))
        logpost2 <- logpost2 + dexp(exp(q[(1 + 3 * ntypes)]), rate = priorcontrol$deltaprior, log = TRUE) + 
            log(exp(q[(1 + 3 * ntypes)]))
        
        ProposedH <- -logpost2 + 0.5 * sum((p^2)/(mass_diag))
        
        ap <- -ProposedH + CurrentH
        alpha <- min(1, exp(ap))
        u <- runif(1)
        if (is.na(alpha)) {
            alpha <- 0
        }
        if (u < alpha) {
            output <- list(logcuralpha = q[1:(1 * ntypes)], logcurbeta = q[(1 * ntypes + 1):(2 * 
                ntypes)], logcurmu = q[(2 * ntypes + 1):(3 * ntypes)], logcurdelta = q[(1 + 3 * 
                ntypes)], logcurgamma = NULL)
            return(output)
        } else {
            output <- list(logcuralpha = logcuralpha, logcurbeta = logcurbeta, logcurmu = logcurmu, 
                logcurdelta = logcurdelta, logcurgamma = NULL)
            return(output)
        }
    }
}


