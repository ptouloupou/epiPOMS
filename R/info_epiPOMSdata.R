#' Modify and extract information from partially observed multi-strain epidemic
#' data.
#'
#' @description
#' This function modifies and extracts general information about an object
#' created by the \code{\link{as_epiPOMSdata}} function.  The output of this
#' function provides the right format and information required to be used in
#' other functions in the package.
#'
#' @param
#' epidata An object of class \sQuote{epiPOMSdata}, produced from the
#' \code{\link{as_epiPOMSdata}} function.
#'
#' @return
#' A list is returned with the following components:
#' \describe{
#' \item{ntypes}{Total number of different strains in the study,
#' \eqn{n_g}.}
#' \item{obserdata}{A list of the observed test results similar to
#' \code{epidata}, but with "-" replaced by 0 and "+" by \eqn{(n_g + 1)}.}
#' \item{ngroups}{Total number of groups in the study, \eqn{P}.}
#' \item{ninds}{A vector of length \eqn{P}, with each
#' element \eqn{p} corresponding to the number of individuals in
#' group \eqn{p}.}
#' \item{tmax}{A vector of length \eqn{P}, with each
#' element \eqn{p} corresponding to the time point that the last sample
#' was collected from the individuals in group \eqn{p}.}
#' \item{indNA}{A vector of length \eqn{P}, with each element \eqn{p}
#' corresponding to the index of the individual that withdrawn from
#' the study before its completion in group \eqn{p}, if any.
#' Otherwise a value of 0 is given.}
#' \item{tmaxNA}{A list containing \eqn{P} vectors of length equal to
#' \code{ninds}\eqn{[p]}, for each \eqn{p = 1, 2, \ldots, P}. Each
#' element of the vector contains the time point that the last sample
#' was collected from the corresponding individual in the group, i.e.
#' if an individual withdrawn from the study before its completion then
#' this value is less than the time of the last observation in the group.}
#' }
#'
#' @details
#' Only works for objects of class \sQuote{epiPOMSdata}. This function is
#' used in conjunction with \code{\link{epiPOMS_mcmc}},
#' \code{\link{plot.epiPOMSmcmc}} and \code{\link{mcmc_initpars}} functions.
#'
#' Note that, in the event where an individual withdrawn from the study, the
#' time of the individual's drop-out is recorded as its last non-missing
#' group pre-specified observation time.
#'
#' @author
#' Panayiota Touloupou
#'
#' @seealso
#' \code{\link{as_epiPOMSdata}} for producing an object of class
#' \sQuote{epiPOMSdata}, \code{\link{epiPOMS_mcmc}} for performing
#' inference in partially observed multi-strain epidemic models and
#' \code{\link{plot.epiPOMSmcmc}} for plotting the outputs of an
#' \sQuote{epiPOMSmcmc} object.
#'
#' @examples
#' # Load E. coli O157:H7 data
#' data(Ecoliepidata)
#' # Load your own dataframe using myData <- read.csv('myDataFile.csv').
#'
#' # Format the data
#' epidata <- as_epiPOMSdata(Ecoliepidata)
#'
#' # Modify and extract information
#' infoepidata <- info_epiPOMSdata(epidata)
#' names(infoepidata)
#'
#' @export
info_epiPOMSdata <- function(epidata) {
    
    if (is(epidata, "epiPOMSdata")) {
        
        DataTest1 <- epidata$epidatatest1
        DataTest2 <- epidata$epidatatest2
        
        ninds <- unlist(lapply(DataTest1, function(x) dim(x)[1]))
        tmax <- unlist(lapply(DataTest1, function(x) dim(x)[2]))
        ngroups <- length(DataTest1)
        
        ntypes1 <- suppressWarnings(as.numeric(unlist(lapply(DataTest1, function(x) max(x, na.rm = T)))))
        ntypes2 <- suppressWarnings(as.numeric(unlist(lapply(DataTest2, function(x) max(x, na.rm = T)))))
        ntypes <- max(c(ntypes1, ntypes2), na.rm = T)
        
        obserdatatest1 <- list()
        obserdatatest2 <- list()
        
        for (p in 1:ngroups) {
            obserdatatest1_group <- matrix(NA, nrow = ninds[p], ncol = tmax[p])
            obserdatatest2_group <- matrix(NA, nrow = ninds[p], ncol = tmax[p])
            for (i in 1:ninds[p]) {
                for (m in 1:tmax[p]) {
                  if (!is.na(DataTest1[[p]][i, m]) & DataTest1[[p]][i, m] == "+") {
                    obserdatatest1_group[i, m] <- as.numeric(ntypes + 1)
                  } else if (!is.na(DataTest1[[p]][i, m]) & DataTest1[[p]][i, m] == "-") {
                    obserdatatest1_group[i, m] <- as.numeric(0)
                  } else {
                    obserdatatest1_group[i, m] <- as.numeric(DataTest1[[p]][i, m])
                  }
                  
                  if (!is.na(DataTest2[[p]][i, m]) & DataTest2[[p]][i, m] == "+") {
                    obserdatatest2_group[i, m] <- as.numeric(ntypes + 1)
                  } else if (!is.na(DataTest2[[p]][i, m]) & DataTest2[[p]][i, m] == "-") {
                    obserdatatest2_group[i, m] <- as.numeric(0)
                  } else {
                    obserdatatest2_group[i, m] <- as.numeric(DataTest2[[p]][i, m])
                  }
                }
            }
            obserdatatest1[[p]] <- obserdatatest1_group
            obserdatatest2[[p]] <- obserdatatest2_group
        }
        
        obsertimes <- list()
        for (p in 1:ngroups) {
            obsertimes_group <- which(!is.na(obserdatatest1[[p]][1, ]))
            for (i in 2:ninds[p]) {
                obsertimes_group <- c(obsertimes_group, which(!is.na(obserdatatest1[[p]][i, ])))
            }
            obsertimes[[p]] <- sort(unique(obsertimes_group))
        }
        
        tmaxNAhelp <- list()
        for (p in 1:ngroups) {
            tmaxNAhelp[[p]] <- suppressWarnings(apply(obserdatatest1[[p]][, obsertimes[[p]]], 1, 
                function(x) min(which(is.na(x))))) - 1
        }
        
        indNA <- rep(NA, ngroups)
        tmaxNA <- list()
        for (p in 1:ngroups) {
            for (i in 1:ninds[p]) {
                if (tmaxNAhelp[[p]][i] == Inf) {
                  tmaxNAhelp[[p]][i] <- tmax[p]
                } else {
                  tmaxNAhelp[[p]][i] <- obsertimes[[p]][tmaxNAhelp[[p]][i]]
                }
            }
            tmaxNA[[p]] <- tmaxNAhelp[[p]]
            indNA[p] <- ifelse(length(which(tmaxNA[[p]] < tmax[p])) > 0, which(tmaxNA[[p]] < tmax[p]), 
                0)
        }
        
    } else {
        stop("info_epiPOMSdata: The epidata must be a class of 'epiPOMSdata'.", call. = FALSE)
    }
    
    list(ntypes = ntypes, obserdata = list(test1 = obserdatatest1, test2 = obserdatatest2), ngroups = ngroups, 
        ninds = ninds, tmax = tmax, indNA = indNA, tmaxNA = tmaxNA)
}

