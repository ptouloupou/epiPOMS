#' Simulate the hidden process (the true carriage states).
#'
#' @description
#' This function simulates the (true) hidden carriage states of each
#' individual per group during its sampling period.
#'
#' @param
#' ntypes An integer greater than 1 providing the number of different strains in
#' the study, \eqn{n_g}.
#' @param
#' ngroups A positive integer providing the total number of groups in the study,
#' \eqn{P}.
#' @param
#' ninds A vector of integers greater than 1 with length \eqn{P}, where
#' element \eqn{p} gives the number of individuals in group \eqn{p}.
#' @param
#' tmax A vector of integers greater than 1 with length \eqn{P}, where
#' element \eqn{p} gives the time point that the last sample was collected from
#' individuals in group \eqn{p}.
#' @param
#' hidpars A list of parameters with which the simulations of the hidden
#' process are to be performed. The names of the components must be "alpha",
#' "beta", "mu", "delta", "gamma" and "nu" (see example). Also, all rates must
#' be specified in days. These components are:
#' \describe{
#' \item{alpha}{A vector of non-negative real numbers, one for each strain,
#' containing the \eqn{\alpha} parameters, i.e. the strain-specific
#' external colonisation rates. This is a vector of \eqn{n_g} values.
#' The first element corresponds to the first strain, the second element to
#' the second strain and so on.}
#' \item{beta}{A vector of non-negative real numbers, one for each strain,
#' containing the \eqn{\beta} parameters, i.e. the strain-specific
#' within-group colonisation rates. This is a vector of \eqn{n_g} values.
#' The first element corresponds to the first strain, the second element to
#' the second strain and so on.}
#' \item{mu}{A vector of non-negative real numbers, one for each strain,
#' containing the \eqn{\mu} parameters, i.e. the strain-specific clearance
#' rates. This is a vector of \eqn{n_g} values.
#' The first element corresponds to the first strain, the second element to
#' the second strain and so on.}
#' \item{delta}{A non-negative real number corresponding to the \eqn{\delta}
#' parameter, i.e. the relative colonisation rate in a carrier versus
#' non-carrier individual.}
#' \item{gamma}{A non-negative real number corresponding to the \eqn{\gamma}
#' parameter, i.e. the relative colonisation rate in smaller versus bigger
#' groups in terms of area (in square meters). If
#' \code{smallgroups = NULL}, i.e there is no difference between groups, set
#' \code{gamma = NULL}.}
#' \item{nu}{A vector of numbers between 0 and 1, containing the \eqn{\nu}
#' parameters, i.e. the probabilities of carriage at the beginning of the study.
#' This is a vector of \eqn{(n_g + 1)} values. The first element corresponds to
#' the non-carriage state. The second and subsequent elements correspond to the
#' carriage of one of the \eqn{n_g} strains. These probabilities should sum up to
#' one.}
#' }
#' @param
#' smallgroups A vector of positive integers, with length equal to the number
#' of groups that are smaller in terms of area (in square meters), containing
#' the indices of the small groups. These values must be between 1 and \eqn{P}.
#' Set \code{smallgroups = NULL} if there is no difference between groups.
#' @param
#' indNA A vector of non-negative integers with length \eqn{P}, where element
#' \eqn{p} gives the index of the individual that withdrawn from the study
#' before its completion in group \eqn{p}, if any. Otherwise a value
#' of 0 should be given. Note that at most one individual from each group can
#' drop-out before the study period ends.
#' @param
#' tNA A vector of positive integers with length \eqn{P}, where element \eqn{p}
#' provides the individual's drop-out time point in group \eqn{p}, if any,
#' otherwise the time of the last observation in group \eqn{p} should be given,
#' \code{tmax}\eqn{[p]}. Note that, if a value is less than the time
#' of the last observation the corresponding value in \code{indNA} must be
#' greater than 0, i.e. if \code{tNA}\eqn{[p]} < \code{tmax}\eqn{[p]}
#' then \code{indNA}\eqn{[p]} > 0.
#'
#' @details
#' Our approach to analyse infection dynamics in longitudinal studies where
#' the process is only partially observed (in the sense that the data are
#' subject to potential testing error due to poor sensitivity of the
#' diagnostic and strain procedure used) involves assuming that the
#' classifications at an observation time are imperfect measures of an underlying
#' true (hidden) epidemic process. Moreover, our approach accounts for missing
#' observations, for example, sparse sampling intervals or individual dropouts.
#'
#' The \code{\link{hidstate_sim}} function simulates the true unobserved
#' (hidden) disease process within each group (the carriage process is
#' assumed to be independent across groups), which is modelled as a multi-state,
#' discrete time non homogeneous Markov model. The possible states include being
#' a non-carrier (state 0) or being a carrier of one of the \eqn{n_g}
#' strains (states \eqn{1, 2, \ldots, n_g}), where 1 to
#' \eqn{(n_g - 1)} refers to carriage of one of the main strains
#' and \eqn{n_g} to carriage of the remaining strains (pooled group).
#' A detailed description of this model is found in
#' \insertCite{Touloupou2020;textual}{epiPOMS}.
#'
#' @return
#' A list containing \eqn{P} matrices (one for each group), with
#' \code{ninds}\eqn{[p]} rows and \code{tmax}\eqn{[p]} columns
#' for each \eqn{p = 1, 2, \ldots, P}. More specifically, if \eqn{h} is the
#' returned list, \eqn{h[[p]]} is the matrix of group \eqn{p} containing the
#' simulated carriage states (of the unobserved process) of each individual in
#' the group at each time point over the sampling period, where the time interval
#' is equal to one day. Therefore, \eqn{h[[p]][c, t]} provides the carriage state
#' of individual \eqn{c = 1, 2, \ldots,} \code{ninds}\eqn{[p]} at time
#' \eqn{t = 1, 2, \ldots,} \code{tmax}\eqn{[p]}, with \code{NA} entries
#' corresponding to missing values, for example, if the individual withdrawn from
#' the study before its completion.
#'
#' @author
#' Panayiota Touloupou, Simon Spencer
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{obserdata_sim}} for simulating partially observed
#' multi-strain epidemic data.
#'
#' @examples
#' # Simulate the hidden carriage state process of individuals housed in
#' # 15 different groups of size varied from 4 to 9. The possible states
#' # include being a non-carrier (state 0) or being a carrier of one of
#' # 5 strains (states 1, 2, 3, 4 and 5).
#' set.seed(1)
#' ntypes <- 5
#' ngroups <- 15
#' ninds <- sample(4:9, size = ngroups, replace = TRUE)
#' tmax <- rep(100, ngroups)
#' hidpars <- list(alpha = c(sample(c(0.0015, 0.001), ntypes - 1, replace = TRUE),
#'      0.002), beta = sample(c(0.005, 0.01), ntypes, replace = TRUE), mu =
#'      sample(c(0.1, 0.15), ntypes, replace = TRUE), delta = 0.5, gamma = 2,
#'      nu = c(0.9, rep((1-0.9)/(ntypes+1), (ntypes-1)), 2*(1-0.9)/(ntypes+1)))
#' smallgroups <- sample(1:ngroups, size = 10)
#' indNA <- c(rep(0, ngroups-1), 1)
#' tNA <- c(rep(100, ngroups-1), 70)
#' hiddenprocess <- hidstate_sim(ntypes = ntypes, ngroups = ngroups,
#'      ninds = ninds, tmax = tmax, hidpars = hidpars, smallgroups =
#'      smallgroups, indNA = indNA, tNA = tNA)
#'
#' @export
hidstate_sim <- function(ntypes, ngroups, ninds, tmax, hidpars, smallgroups, indNA, tNA) {
    
    # Error checks for input arguments
    check <- ntypes
    if (is.null(check)) {
        stop("hidstate_sim: The ntypes has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("hidstate_sim: The ntypes must be an integer greater than one.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("hidstate_sim: The ntypes must be an integer greater than one.", call. = FALSE)
    } else if ((check <= 1) | (check != round(ntypes))) {
        stop("hidstate_sim: The ntypes must be an integer greater than one.", call. = FALSE)
    }
    
    check <- ngroups
    if (is.null(check)) {
        stop("hidstate_sim: The ngroups has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("hidstate_sim: The ngroups must be a positive integer.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("hidstate_sim: The ngroups must be a positive integer.", call. = FALSE)
    } else if ((check <= 0) | (check != round(check))) {
        stop("hidstate_sim: The ngroups must be a positive integer.", call. = FALSE)
    }
    
    check <- ninds
    if (is.null(check)) {
        stop("hidstate_sim: The ninds has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The ninds must be a vector.", call. = FALSE)
    } else if (length(check) != ngroups) {
        stop("hidstate_sim: Length of ninds is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of ninds must be integers greater than one.", call. = FALSE)
    } else if (any(check <= 1) | any(check != round(check))) {
        stop("hidstate_sim: The entries of ninds must be integers greater than one.", call. = FALSE)
    }
    
    check <- tmax
    if (is.null(check)) {
        stop("hidstate_sim: The tmax has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The tmax must be a vector.", call. = FALSE)
    } else if (length(check) != ngroups) {
        stop("hidstate_sim: Length of tmax is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of tmax must be integers greater than one.", call. = FALSE)
    } else if (any(check <= 1) | any(check != round(check))) {
        stop("hidstate_sim: The entries of tmax must be integers greater than one.", call. = FALSE)
    }
    
    check <- hidpars$alpha
    if (is.null(check)) {
        stop("hidstate_sim: The alpha at hidpars has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The alpha at hidpars must be a vector.", call. = FALSE)
    } else if (length(check) != ntypes) {
        stop("hidstate_sim: Length of alpha at hidpars is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of alpha at hidpars must be non-negative reals.", call. = FALSE)
    } else if (any(check < 0)) {
        stop("hidstate_sim: The entries of alpha at hidpars must be non-negative reals.", call. = FALSE)
    }
    
    check <- hidpars$beta
    if (is.null(check)) {
        stop("hidstate_sim: The beta at hidpars has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The beta at hidpars must be a vector.", call. = FALSE)
    } else if (length(check) != ntypes) {
        stop("hidstate_sim: Length of beta at hidpars is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of beta at hidpars must be non-negative reals.", call. = FALSE)
    } else if (any(check < 0)) {
        stop("hidstate_sim: The entries of beta at hidpars must be non-negative reals.", call. = FALSE)
    }
    
    check <- hidpars$mu
    if (is.null(check)) {
        stop("hidstate_sim: The mu at hidpars has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The mu at hidpars must be a vector.", call. = FALSE)
    } else if (length(check) != ntypes) {
        stop("hidstate_sim: Length of mu at hidpars is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of mu at hidpars must be non-negative reals.", call. = FALSE)
    } else if (any(check < 0)) {
        stop("hidstate_sim: The entries of mu at hidpars must be non-negative reals.", call. = FALSE)
    }
    
    check <- hidpars$delta
    if (is.null(check)) {
        stop("hidstate_sim: The delta at hidpars has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("hidstate_sim: The delta at hidpars must be a non-negative real.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("hidstate_sim: The delta at hidpars must be a non-negative real.", call. = FALSE)
    } else if ((check < 0)) {
        stop("hidstate_sim: The delta at hidpars must be a non-negative real.", call. = FALSE)
    }
    
    check <- hidpars$gamma
    if (!is.null(check)) {
        if ((length(check) != 1) | is.list(check)) {
            stop("hidstate_sim: The gamma at hidpars must be a non-negative real.", call. = FALSE)
        } else if (!is.numeric(check)) {
            stop("hidstate_sim: The gamma at hidpars must be a non-negative real.", call. = FALSE)
        } else if ((check < 0)) {
            stop("hidstate_sim: The gamma at hidpars must be a non-negative real.", call. = FALSE)
        }
    }
    
    if ((is.null(smallgroups)) != (is.null(check))) {
        stop("hidstate_sim: The smallgroups must be NULL if gamma at hidpars is NULL, and via versa.", 
            call. = FALSE)
    }
    
    check <- hidpars$nu
    if (is.null(check)) {
        stop("hidstate_sim: The nu at hidpars has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The nu at hidpars must be a vector.", call. = FALSE)
    } else if (length(check) != (ntypes + 1)) {
        stop("hidstate_sim: Length of nu at hidpars is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of nu at hidpars must be between zero and one.", call. = FALSE)
    } else if (any(check < 0) | any(check > 1)) {
        stop("hidstate_sim: The entries of nu at hidpars must be between zero and one.", call. = FALSE)
    } else if (sum(check) > 1) {
        stop("hidstate_sim:  The entries of nu at hidpars must sum to one.", call. = FALSE)
    }
    
    check <- smallgroups
    if (!is.null(check)) {
        if ((!is.vector(check)) | is.list(check)) {
            stop("hidstate_sim: The smallgroups must be a vector or NULL.", call. = FALSE)
        } else if (length(check) > ngroups) {
            stop("hidstate_sim: Length of smallgroups is not compatible.", call. = FALSE)
        } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
            stop("hidstate_sim: The entries of smallgroups must be integers between 1 and ngroups or NULL.", 
                call. = FALSE)
        } else if (any(check < 1) | any(check > ngroups) | any(check != round(check))) {
            stop("hidstate_sim: The entries of smallgroups must be integers between 1 and ngroups or NULL.", 
                call. = FALSE)
        }
    }
    
    check <- indNA
    if (is.null(check)) {
        stop("hidstate_sim: The indNA has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The indNA must be a vector.", call. = FALSE)
    } else if (length(check) != ngroups) {
        stop("hidstate_sim: Length of indNA is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of indNA must be integers between 0 and ninds in the group.", 
            call. = FALSE)
    } else if (any(check < 0) | any(check > ninds) | any(check != round(check))) {
        stop("hidstate_sim: The entries of indNA must be integers between 0 and ninds in the group.", 
            call. = FALSE)
    }
    
    check <- tNA
    if (is.null(check)) {
        stop("hidstate_sim: The tNA has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("hidstate_sim: The tNA must be a vector.", call. = FALSE)
    } else if (length(check) != ngroups) {
        stop("hidstate_sim: Length of tNA is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("hidstate_sim: The entries of tNA must be integers between 1 and tmax of the group.", 
            call. = FALSE)
    } else if (any(check <= 0) | any(check > tmax) | any(check != round(check))) {
        stop("hidstate_sim: The entries of tNA must be integers between 1 and tmax of the group.", 
            call. = FALSE)
    } else if (any(indNA[which(tNA < tmax)] == 0)) {
        stop("hidstate_sim: The entries of indNA in which tNA < tmax must be greater than zero.", 
            call. = FALSE)
    }
    # End of error checks for input arguments
    
    alltypes <- 0
    while (alltypes != (ntypes + 1)) {
        hiddenprocess <- list()
        
        for (i in 1:ngroups) {
            hiddenprocesshelp <- matrix(NA, ninds[i], tmax[i])
            gamma_new <- ifelse(i %in% smallgroups, hidpars$gamma, 1)
            
            hiddenprocesshelp[, 1] <- apply(rmultinom(n = ninds[i], size = 1, prob = hidpars$nu), 
                2, function(x) which(x == 1)) - 1
            
            for (j in 2:tNA[i]) {
                It <- rep(0, ntypes)
                zInfec <- subset(hiddenprocesshelp[, j - 1], hiddenprocesshelp[, j - 1] > 0)
                tI <- table(zInfec)
                thesis <- as.numeric(names(tI))
                It[thesis] <- tI
                Rec <- hidpars$mu
                Infec <- hidpars$alpha + It * hidpars$beta * gamma_new
                
                Trans <- matrix(NA, ntypes + 1, ntypes + 1)
                Trans[1, 1] <- exp(-sum(Infec))
                probhelp <- (1 - exp(-sum(Infec)))/(sum(Infec))
                Trans[1, -1] <- (Infec) * probhelp
                
                for (ts in 1:ntypes) {
                  xtr <- c(1:ntypes)[-ts]
                  probhelp <- (1 - (exp(-hidpars$mu[ts] - hidpars$delta * sum(Infec[xtr]))))/(hidpars$mu[ts] + 
                    hidpars$delta * sum(Infec[xtr]))
                  Trans[ts + 1, 1] <- (hidpars$mu[ts]) * probhelp
                  
                  for (tr in xtr) {
                    Trans[ts + 1, tr + 1] <- (hidpars$delta * Infec[tr]) * probhelp
                  }
                  Trans[ts + 1, ts + 1] <- ifelse((1 - sum(Trans[ts + 1, -(ts + 1)])) < 0, 0, 1 - 
                    sum(Trans[ts + 1, -(ts + 1)]))
                }
                
                for (ind in 1:ninds[i]) {
                  for (ts in 0:ntypes) {
                    
                    if (hiddenprocesshelp[ind, j - 1] == ts) {
                      hiddenprocesshelp[ind, j] <- which(rmultinom(n = 1, size = 1, prob = Trans[ts + 
                        1, ]) == 1) - 1
                    }
                  }
                }
                
            }
            
            if (tNA[i] != tmax[i]) {
                for (j in (tNA[i] + 1):tmax[i]) {
                  It <- rep(0, ntypes)
                  zInfec <- subset(hiddenprocesshelp[-indNA[i], j - 1], hiddenprocesshelp[-indNA[i], 
                    j - 1] > 0)
                  tI <- table(zInfec)
                  thesis <- as.numeric(names(tI))
                  It[thesis] <- tI
                  
                  Rec <- hidpars$mu
                  Infec <- hidpars$alpha + It * hidpars$beta * gamma_new
                  
                  Trans <- matrix(NA, ntypes + 1, ntypes + 1)
                  Trans[1, 1] <- exp(-sum(Infec))
                  probhelp <- (1 - exp(-sum(Infec)))/sum(Infec)
                  Trans[1, -1] <- Infec * probhelp
                  
                  for (ts in 1:ntypes) {
                    xtr <- c(1:ntypes)[-ts]
                    probhelp <- (1 - (exp(-hidpars$mu[ts] - hidpars$delta * sum(Infec[xtr]))))/(hidpars$mu[ts] + 
                      hidpars$delta * sum(Infec[xtr]))
                    Trans[ts + 1, 1] <- (hidpars$mu[ts]) * probhelp
                    
                    for (tr in xtr) {
                      Trans[ts + 1, tr + 1] <- (hidpars$delta * Infec[tr]) * probhelp
                    }
                    Trans[ts + 1, ts + 1] <- ifelse((1 - sum(Trans[ts + 1, -(ts + 1)])) < 0, 0, 
                      1 - sum(Trans[ts + 1, -(ts + 1)]))
                    
                  }
                  
                  
                  for (ind in (1:ninds[i])[-indNA[i]]) {
                    for (ts in 0:ntypes) {
                      
                      if (hiddenprocesshelp[ind, j - 1] == ts) {
                        hiddenprocesshelp[ind, j] <- which(rmultinom(n = 1, size = 1, prob = Trans[ts + 
                          1, ]) == 1) - 1
                      }
                    }
                  }
                  
                }
                
                
            }
            hiddenprocess[[i]] <- hiddenprocesshelp
        }
        alltypes <- length(unique(unlist(lapply(hiddenprocess, function(x) unique(na.omit(as.vector(x)))))))
    }
    return(hiddenprocess)
}

