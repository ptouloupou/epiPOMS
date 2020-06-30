#' Simulate partially observed multi-strain epidemic data.
#'
#' @description
#' This function allows the user to generate simulations of the measurement
#' process (the observed test results) of a population of individuals which is
#' partitioned into groups. This function returns an object of class
#' \sQuote{epiPOMSdata}.
#'
#' @inheritParams hidstate_sim
#' @param
#' obserpars A vector of 5 values, with which the simulations of the observed
#' process are to be performed, containing the following parameters in order:
#' \eqn{\theta_1}, \eqn{\theta_2}, \eqn{\theta_C}, \eqn{\theta_S} and
#' \eqn{\theta_P}, where \eqn{\theta_1} and \eqn{\theta_2} denote the test
#' sensitivities. Given that a test is positive, \eqn{\theta_C} denotes
#' the probability of correctly identifying a common strain, \eqn{\theta_S}
#' is the probability of misclassifying a common strain with a different common
#' strain and \eqn{\theta_P} the probability that a strain of pooled type is
#' classified as a common strain. All parameters must be between 0 and 1 and
#' \eqn{\theta_C + \theta_S \leq 1}{\theta_C + \theta_S \le 1}. See details,
#' below, for more information about test classification.
#' @param
#' obsertimes A list containing \eqn{P} vectors of positive integers
#' giving the observation/sampling times, i.e. the times that individuals
#' were tested in each group \eqn{p = 1, 2, \dots, P}, with entries ranging
#' from one to \code{tmax}\eqn{[p]}. Each time vector must be in ascending
#' order, have length greater than 1 and first element equal to 1 (corresponding
#' to the time that the first sample was collected from the individuals in that
#' group). Moreover, the last element of each vector \eqn{p} must be equal to
#' \code{tmax}\eqn{[p]} (see example).
#' @param
#' npostyped A positive integer providing the number of positive samples
#' (from either test) that will be randomly selected to be strain typed from
#' each group of the study. When less than \code{npostyped} positive samples
#' obtained within a group, all the available samples are chosen to be typed.
#'
#' @details
#' Our approach to analyse infection dynamics in longitudinal studies where
#' the process is only partially observed, in the sense that the data are
#' subject to potential testing error due to poor sensitivity of the
#' diagnostic and strain procedure used, involves assuming that the
#' classifications at an observation time are imperfect measures of an underlying
#' true (hidden) epidemic process. Moreover, our approach accounts for missing
#' observations, for example, sparse sampling intervals or individual dropouts.
#'
#' As noted above, the underlying carriage process is not directly observed.
#' Instead, for each individual we obtain the results of 2 diagnostic tests
#' that are taken on pre-specified discrete times in each group. Function
#' \code{\link{obserdata_sim}} simulates the observed process within
#' each group for each time point over the sampling period. It should be noted
#' that the interval between two time points is one day. At each sampling time
#' the observed states are generated conditional on the true disease
#' states, as obtained by function \code{\link{hidstate_sim}}, according
#' to a misspecification matrix (one for each test). Once the complete data set
#' of fully strain typed test results is generated, we randomly choose
#' \code{npostyped} strain typed observations per group to remain typed and set
#' the remaining as untyped positives. When less than \code{npostyped} typed
#' samples occured within a group, all samples remain typed.
#'
#' To sum up, at an observation time a test result is classified as "-", if
#' it is negative, "+" if it is positive but not chosen for strain typing,
#' otherwise, strain number is given ranging from 1 to \eqn{n_g}, where
#' 1 to \eqn{(n_g - 1)} refers to carriage of one of the main strains and
#' \eqn{n_g} to carriage of the remaining strains (pooled group). \code{NA}
#' entries corresponds to missing values, either because the individual is
#' not tested, for example due to sparse sampling intervals, or withdraw from
#' the study before its completion. A detailed description is found in
#' \insertCite{Touloupou2020;textual}{epiPOMS}.
#'
#' @return
#' An object of class \sQuote{epiPOMSdata} is returned containing the simulated
#' observed test results of all individuals. For a precise description,
#' see \code{\link{as_epiPOMSdata}}.
#'
#' @author
#' Panayiota Touloupou, Simon Spencer
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for performing inference in partially
#' observed multi-strain epidemic models.
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
#' # Plot the simulated epidemic data
#' plot(epidata)
#'
#' @export
obserdata_sim <- function(ntypes, ngroups, ninds, tmax, hidpars, obserpars, smallgroups, obsertimes, 
    indNA, tNA, npostyped) {
    
    # Error checks for input arguments
    check <- obserpars
    if (is.null(check)) {
        stop("obserdata_sim: The obserpars has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("obserdata_sim: The obserpars must be a vector.", call. = FALSE)
    } else if (length(check) != 5) {
        stop("obserdata_sim: Length of obserpars is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("obserdata_sim: The entries of obserpars must be between zero and one.", call. = FALSE)
    } else if (any(check < 0) | any(check > 1)) {
        stop("obserdata_sim: The entries of obserpars must be between zero and one.", call. = FALSE)
    } else if (sum(check[3:4]) > 1) {
        stop("obserdata_sim: The sum of third and forth elements of obserpars must be less than one.", 
            call. = FALSE)
    }
    
    check <- obsertimes
    if (is.null(check)) {
        stop("obserdata_sim: The obsertimes has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("obserdata_sim: The obsertimes must be a list.", call. = FALSE)
    } else if (length(check) != ngroups) {
        stop("obserdata_sim: Length of obsertimes is not compatible", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) ((!is.vector(x)) | is.list(x)))))) {
        stop("obserdata_sim: The obsertimes must contain P vectors of length greater than one.", 
            call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (length(x) == 1))))) {
        stop("obserdata_sim: The obsertimes must contain P vectors of length greater than one.", 
            call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) suppressWarnings(any(is.na(as.numeric(x)))))))) {
        stop("obserdata_sim: The entries of obsertimes must be positive integers.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (any(x <= 0) | any(x != round(x))))))) {
        stop("obserdata_sim: The entries of obsertimes must be positive integers.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (sort(x) != x))))) {
        stop("obserdata_sim: The obsertimes must contain P vectors in ascending order.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (x[1] != 1))))) {
        stop("obserdata_sim: The obsertimes must contain P vectors with first element equal to 1.", 
            call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) max(x))) != tmax)) {
        stop("obserdata_sim: The obsertimes must contain P vectors with last element equal to tmax[p].", 
            call. = FALSE)
    }
    
    check <- npostyped
    if (is.null(check)) {
        stop("obserdata_sim: The npostyped has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("obserdata_sim: The npostyped must be a positive integer.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("obserdata_sim: The npostyped must be a positive integer.", call. = FALSE)
    } else if ((check <= 0) | (check != round(check))) {
        stop("obserdata_sim: The npostyped must be a positive integer.", call. = FALSE)
    }
    
    # End of error checks for input arguments
    
    
    alltypes <- 0
    while (alltypes != ntypes) {
        
        hiddenprocess <- hidstate_sim(ntypes = ntypes, ngroups = ngroups, ninds = ninds, tmax = tmax, 
            hidpars = hidpars, smallgroups = smallgroups, indNA = indNA, tNA = tNA)
        
        classmatrices <- classprob(ntypes = ntypes, obserpars = obserpars)
        
        
        Test1 <- list()
        Test2 <- list()
        
        for (p in 1:ngroups) {
            
            Test1_Group <- matrix(NA, nrow = ninds[p], ncol = tmax[p])
            Test2_Group <- matrix(NA, nrow = ninds[p], ncol = tmax[p])
            
            for (i in 1:ninds[p]) {
                for (m in 1:tmax[p]) {
                  for (ts in 0:ntypes) {
                    if (!is.na(hiddenprocess[[p]][i, m]) & (hiddenprocess[[p]][i, m] == ts) & (m %in% 
                      obsertimes[[p]])) {
                      Test1_Group[i, m] <- which(rmultinom(n = 1, size = 1, prob = classmatrices$classprobtest1[-(ntypes + 
                        2), ts + 1]) == 1) - 1
                      Test2_Group[i, m] <- which(rmultinom(n = 1, size = 1, prob = classmatrices$classprobtest2[-(ntypes + 
                        2), ts + 1]) == 1) - 1
                    }
                  }
                }
            }
            
            pos1 <- which(Test1_Group > 0, arr.ind = T)
            pos2 <- which(Test2_Group > 0, arr.ind = T)
            
            if ((nrow(pos1) > 0) & (nrow(pos2) > 0)) {
                
                Test1Help <- Test1_Group
                Test2Help <- Test2_Group
                Test1Help[pos1] <- ntypes + 1
                Test2Help[pos2] <- ntypes + 1
                
                pos12 <- rbind(cbind(pos1, 1), cbind(pos2, 2))
                if (nrow(pos12) < npostyped) {
                  rand <- 1:nrow(pos12)
                } else {
                  rand <- sample(1:dim(pos12)[1], npostyped, replace = F)
                }
                randstrain <- pos12[rand, ]
                
                for (j in 1:nrow(randstrain)) {
                  if (randstrain[j, 3] == 1) {
                    Test1Help[randstrain[j, 1], randstrain[j, 2]] <- Test1_Group[randstrain[j, 
                      1], randstrain[j, 2]]
                  } else {
                    Test2Help[randstrain[j, 1], randstrain[j, 2]] <- Test2_Group[randstrain[j, 
                      1], randstrain[j, 2]]
                  }
                  
                }
                Test1Help[which(Test1Help == 0)] <- "-"
                Test2Help[which(Test2Help == 0)] <- "-"
                Test1Help[which(Test1Help == (ntypes + 1))] <- "+"
                Test2Help[which(Test2Help == (ntypes + 1))] <- "+"
                
                Test1[[p]] <- Test1Help
                Test2[[p]] <- Test2Help
            } else {
                Test1_Group[which(Test1_Group == 0)] <- "-"
                Test2_Group[which(Test2_Group == 0)] <- "-"
                
                Test1[[p]] <- Test1_Group
                Test2[[p]] <- Test2_Group
                
            }
        }
        
        
        alltypes1 <- unique(unlist(lapply(Test1, function(x) unique(na.omit(as.vector(x))))))
        alltypes2 <- unique(unlist(lapply(Test2, function(x) unique(na.omit(as.vector(x))))))
        
        alltypes1 <- alltypes1[which(!(alltypes1 %in% c("-", "+")))]
        alltypes2 <- alltypes2[which(!(alltypes2 %in% c("-", "+")))]
        
        alltypes <- length(unique(c(alltypes1, alltypes2)))
    }
    output <- list(epidatatest1 = Test1, epidatatest2 = Test2)
    class(output) <- "epiPOMSdata"
    return(output)
}

