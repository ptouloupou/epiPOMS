#' Produce test classification matrices.
#'
#' @description
#' Provides the test classification matrices containing the probabilities of
#' correctly or incorrectly classifying the observed outcome given the true
#' carriage state. The output of this function provides the right format
#' required to be used into other functions in the package.
#'
#' @param
#' ntypes An integer greater than 1 providing the number of different strains in
#' the study.
#' @param
#' obserpars A vector of 5 values, containing the following observation model
#' parameters in order: \eqn{\theta_1}, \eqn{\theta_2}, \eqn{\theta_C},
#' \eqn{\theta_S} and \eqn{\theta_P}, where \eqn{\theta_1} and \eqn{\theta_2}
#' denote the test sensitivities. Given that a test is positive, \eqn{\theta_C}
#' denotes the probability of correctly identifying a common strain,
#' \eqn{\theta_S} is the probability of misclassifying a common strain with a
#' different common strain and \eqn{\theta_P} the probability that a strain of
#' pooled type is classified as a common strain. All parameters must be between
#' 0 and 1 and \eqn{\theta_C + \theta_S \leq 1}{\theta_C + \theta_S \le 1}.
#'
#' @details
#' Auxiliary function that can be used to provide the test classification
#' probabilities used by \code{\link{epiPOMS_mcmc}} for performing Bayesian
#' inference in partially observed multi-strain epidemic models.
#'
#' @return
#' A list is returned with the following components:
#' \describe{
#' \item{classprobtest1}{A matrix, with \code{(ntypes + 2)} rows and
#' \code{(ntypes + 1)} columns, containing the classification probabilities of
#' the first test.}
#' \item{classprobtest2}{A matrix with the same structure as
#' \code{classprobtest1}, but for the second test.}
#' }
#'
#' @author
#' Panayiota Touloupou
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for performing inference in partially
#' observed multi-strain epidemic models.
#'
#' @export
classprob <- function(ntypes, obserpars) {
    
    classtest1 <- matrix(0, ntypes + 2, ntypes + 1)
    classtest1[, 1] <- c(1, rep(0, ntypes + 1))
    classtest1[ntypes + 2, -1] <- obserpars[1]
    classtest1[1, -1] <- 1 - obserpars[1]
    classtest1[-c(1, ntypes + 1, ntypes + 2), -c(1, ntypes + 1)] <- obserpars[1] * obserpars[4]/(ntypes - 
        2)
    classtest1[-c(1, ntypes + 1, ntypes + 2), ntypes + 1] <- obserpars[1] * obserpars[5]/(ntypes - 
        1)
    for (ij in 2:ntypes) {
        classtest1[ij, ij] <- obserpars[1] * obserpars[3]
    }
    classtest1[ntypes + 1, -c(1, ntypes + 1)] <- obserpars[1] * (1 - obserpars[3] - obserpars[4])
    classtest1[ntypes + 1, ntypes + 1] <- obserpars[1] * (1 - obserpars[5])
    
    classtest2 <- matrix(0, ntypes + 2, ntypes + 1)
    classtest2[, 1] <- c(1, rep(0, ntypes + 1))
    classtest2[ntypes + 2, -1] <- obserpars[2]
    classtest2[1, -1] <- 1 - obserpars[2]
    classtest2[-c(1, ntypes + 1, ntypes + 2), -c(1, ntypes + 1)] <- obserpars[2] * obserpars[4]/(ntypes - 
        2)
    classtest2[-c(1, ntypes + 1, ntypes + 2), ntypes + 1] <- obserpars[2] * obserpars[5]/(ntypes - 
        1)
    for (ij in 2:ntypes) {
        classtest2[ij, ij] <- obserpars[2] * obserpars[3]
    }
    classtest2[ntypes + 1, -c(1, ntypes + 1)] <- obserpars[2] * (1 - obserpars[3] - obserpars[4])
    classtest2[ntypes + 1, ntypes + 1] <- obserpars[2] * (1 - obserpars[5])
    
    
    output <- list(classprobtest1 = classtest1, classprobtest2 = classtest2)
    return(output)
}

