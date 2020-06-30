#' Convert multi-strain epidemic dataframe into an \sQuote{epiPOMSdata} object.
#'
#' @description
#' This function allows the user to generate objects of class
#' \sQuote{epiPOMSdata}. The output of this function provides the right format
#' for input into other functions in the package.
#'
#' @param
#' epidata A dataframe consisting of multi-strain epidemic data obtained
#' from a longitudinal study of a population of individuals which is divided
#' into groups. For each individual, the results of two diagnostic tests taken
#' on pre-specified discrete times are provided. More specifically:
#' \describe{
#' \item{Column 1}{Gives the time at which each observation is made: an integer
#' ranging from one to \eqn{T}, where \eqn{T > 1} is the maximum observation time
#' point of the study. See details, below, for more information.}
#' \item{Column 2}{Gives the individual index from which each observation is made:
#' an integer starting from one to \eqn{C}, where \eqn{C > 1} is the maximum
#' number of individuals per group in the study. See details, below, for more
#' information.}
#' \item{Column 3}{Gives the group index of the individual from which each
#' observation is made: an integer starting from one to \eqn{P}, where \eqn{P} is
#' the total number of groups in the study.}
#' \item{Columns 4 and 5}{Give the observed result of the first and second
#' diagnostic test, respectively, of each observation: "-", if the test result is
#' negative, "+" if it is positive but not chosen for strain typing, otherwise,
#' strain number is given (starting from one to \eqn{n_g > 1}, where
#' \eqn{n_g} is the total number of different strains in the study). See data
#' restrictions in details, below, for more information.}
#' }
#'
#' @details
#' The entries of column 1 must start from one for each individual in a group,
#' with one being the time that the first sample was collected from the
#' individuals in that particular group, i.e. initial time. Moreover, individuals
#' within a group must be tested more than once, at any later time than the
#' initial. The entries of column 2 must start from one up to the total number of
#' individuals in each group (which must be greater than 1), including all
#' positive numbers between the specific range. The same restrictions apply to
#' group and strain indices.
#'
#' See \code{\link{Ecoliepidata}} as an example of such a dataframe.
#'
#' \strong{Data restrictions:}
#'
#' The data should contain observations from a population of individuals
#' partitioned into groups (possibly of various sizes) with more than one
#' members, so that each individual belongs to only one group for the entire
#' period of the study. For example, these groups can be households within a
#' community or pens within a feedlot. Individuals within the same group are
#' tested at common time points (individuals in other groups may be tested on
#' different times) and on these observation times we obtain the results of two
#' diagnostic tests. Nevertheless, individuals are allowed to miss a
#' pre-specified observation time or withdraw from the study or no longer
#' be able to followed up for some other reason. However, we assume that at most
#' one individual from each group can drop-out before the study period ends.
#' The data provided by the user must not include any missing observations, but
#' only the times that samples are collected from each individual.
#'
#' Finally, at each sampling time a test result is classified as "-", if
#' it is negative, "+" if it is positive but not chosen for strain typing,
#' otherwise, strain number is given ranging from 1 to \eqn{n_g > 1}, where
#' 1 to \eqn{(n_g - 1)} refers to carriage of one of the
#' main strains and \eqn{n_g} to carriage of the remaining strains.
#' This formulation implies that the remaining strains are treated as a single
#' group, referred as the "pooled" group, and assumed to be of the same type
#' \eqn{n_g}.
#'
#' @return
#' An object of class \sQuote{epiPOMSdata} is returned containing the following:
#' \describe{
#' \item{epidatatest1}{A list containing \eqn{P} matrices with the observed
#' results of the first test, one for each group \eqn{p = 1, 2, \ldots, P}.
#' Each matrix has \eqn{C^{[p]}}{C^[p]} rows and \eqn{T^{[p]}}{T^[p]}
#' columns, where
#' \eqn{C^{[p]} \subseteq \{1, 2, \ldots, C\}}{C^[p] ⊆ \{1, 2, \ldots,
#' C\}} and \eqn{T^{[p]} \subseteq \{1, 2, \ldots, T\}}{T^[p] ⊆ \{1, 2,
#' \ldots, T\}} denote the total number of individuals and the time that the
#' last sample was collected in group \eqn{p}, respectively. More specifically,
#' if \eqn{d} is the returned list, \eqn{d[[p]]} is the matrix of group \eqn{p}
#' containing the test result of each individual in the group for each time
#' point, where the time interval is equal to one day. Therefore,
#' \eqn{d[[p]][c, t]} provides the observation of individual
#' \eqn{c = 1, 2, \ldots, C^{[p]}}{c = 1, 2, \ldots, C^[p]} at time
#' \eqn{t = 1, 2, \ldots, T^{[p]}}{t = 1, 2, \ldots, T^[p]}, with \code{NA}
#' entries corresponding to missing observations that may have occurred due to,
#' for example, sparse sampling intervals or individual dropouts.}
#' \item{epidatatest2}{A list with the same structure as \code{epidatatest1},
#' but for the second test.}
#' }
#'
#' @author
#' Panayiota Touloupou
#'
#' @seealso
#' \code{\link{plot.epiPOMSdata}} for plotting an object of class
#' \sQuote{epiPOMSdata}, \code{\link{epiPOMS_mcmc}} for performing
#' inference in partially observed multi-strain epidemic models.
#'
#' @examples
#' # Load E. coli O157:H7 data
#' data(Ecoliepidata)
#' # Load your own dataframe using myData <- read.csv('myDataFile.csv').
#' # See data restrictions.
#'
#' # Format the data
#' epiPOMSdata <- as_epiPOMSdata(Ecoliepidata)
#' names(epiPOMSdata)
#'
#' # Plot the E. coli O157:H7 data
#' plot(epiPOMSdata)
#'
#' @export
as_epiPOMSdata <- function(epidata) {
    
    # Error checks for input arguments
    if (is.null(epidata)) {
        stop("as_epiPOMSdata: The epidata has to be specified.", call. = FALSE)
    } else if (!is.data.frame(epidata)) {
        stop("as_epiPOMSdata: The epidata has to be specified as a dataframe.", call. = FALSE)
    } else if (ncol(epidata) != 5) {
        stop("as_epiPOMSdata: The epidata must have 5 columns.", call. = FALSE)
    }
    
    check1 <- epidata[which(!(epidata[, 4] %in% c("-", "+"))), 4]
    check1 <- suppressWarnings(as.numeric(levels(check1)[check1]))
    check2 <- epidata[which(!(epidata[, 5] %in% c("-", "+"))), 5]
    check2 <- suppressWarnings(as.numeric(levels(check2)[check2]))
    check <- c(check1, check2)
    
    if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("as_epiPOMSdata: The entries of columns 4 and 5 of epidata must be \"+\" or \"-\" or positive integers.", 
            call. = FALSE)
    } else if (any(check <= 0) | any(check != round(check))) {
        stop("as_epiPOMSdata: The entries of columns 4 and 5 of epidata must be \"+\" or \"-\" or positive integers.", 
            call. = FALSE)
    } else if (length(check) == 0) {
        stop("as_epiPOMSdata: The entries of columns 4 and 5 of epidata must also contain positive integers and not only \"+\" and \"-\".", 
            call. = FALSE)
    } else if (max(check) == 1) {
        stop("as_epiPOMSdata: Number of different strains must be greater than one.", call. = FALSE)
    } else if (length(unique(check)) != max(check)) {
        stop("as_epiPOMSdata: Strains must be indexed from one to the total number of strains in the study.", 
            call. = FALSE)
    } else if (any(sort(unique(check)) != c(1:max(check)))) {
        stop("as_epiPOMSdata: Groups must be indexed from one to the total number of groups.", 
            call. = FALSE)
    }
    
    check <- epidata[, 1:3]
    if (any(apply(check, MARGIN = 2, FUN = function(x) suppressWarnings(any(is.na(as.numeric(x))))))) {
        stop("as_epiPOMSdata: The entries of columns 1,2 and 3 of epidata must be positive integers.", 
            call. = FALSE)
    } else if (any(check <= 0) | any(check != round(check))) {
        stop("as_epiPOMSdata: The entries of columns 1,2 and 3 of epidata must be positive integers.", 
            call. = FALSE)
    }
    
    ngroups <- max(epidata[, 3])
    
    if (length(unique(epidata[, 3])) != ngroups) {
        stop("as_epiPOMSdata: Groups must be indexed from one to the total number of groups.", 
            call. = FALSE)
    } else if (any(sort(unique(epidata[, 3])) != c(1:ngroups))) {
        stop("as_epiPOMSdata: Groups must be indexed from one to the total number of groups.", 
            call. = FALSE)
    }
    
    for (i in 1:ngroups) {
        epigroup <- subset(epidata, epidata[, 3] == i)
        maxtime <- max(epigroup[, 1])
        maxind <- max(epigroup[, 2])
        
        if (maxtime == 1) {
            stop("as_epiPOMSdata: The individuals in a group must be sampled at least two times.", 
                call. = FALSE)
        }
        
        if (length(unique(epigroup[, 2])) != maxind) {
            stop("as_epiPOMSdata: Individuals per group must be indexed from one to the total number of individuals in that group.", 
                call. = FALSE)
        } else if (any(sort(unique(epigroup[, 2])) != c(1:maxind))) {
            stop("as_epiPOMSdata: Individuals per group must be indexed from one to the total number of individuals in that group.", 
                call. = FALSE)
        }
        
        if (maxind == 1) {
            stop("as_epiPOMSdata: The groups must have more than one members.", call. = FALSE)
        }
        
        if (sum(maxtime == epigroup[, 1]) < (maxind - 1)) {
            stop("as_epiPOMSdata: Only one individual from each group can withdraw from the study.", 
                call. = FALSE)
        }
    }
    
    # End of error checks for input arguments
    
    colnames(epidata) <- c("Time", "Individual", "Group", "ResultTest1", "ResultTest2")
    
    ngroups <- max(epidata$Group)
    Test1 <- list()
    Test2 <- list()
    for (i in 1:ngroups) {
        ninds <- max(epidata$Individual[(epidata$Group == i)])
        tmax <- max(epidata$Time[(epidata$Group == i)])
        MatrixHelpTest1 <- matrix(NA, nrow = ninds, ncol = tmax)
        MatrixHelpTest2 <- matrix(NA, nrow = ninds, ncol = tmax)
        for (j in 1:ninds) {
            tmaxIndiv <- epidata$Time[(epidata$Group == i) & (epidata$Individual == j)]
            TestIndiv <- as.character(epidata$ResultTest1[(epidata$Group == i) & (epidata$Individual == 
                j)])
            MatrixHelpTest1[j, tmaxIndiv] <- TestIndiv
            TestIndiv <- as.character(epidata$ResultTest2[(epidata$Group == i) & (epidata$Individual == 
                j)])
            MatrixHelpTest2[j, tmaxIndiv] <- TestIndiv
        }
        
        Test1[[i]] <- MatrixHelpTest1
        Test2[[i]] <- MatrixHelpTest2
    }
    
    output <- list(epidatatest1 = Test1, epidatatest2 = Test2)
    class(output) <- "epiPOMSdata"
    return(output)
}




