#' Plot partially observed multi-strain epidemic data.
#'
#' @description
#' This function plots multi-strain epidemic longitudinal data consists
#' of a set of result sequences from two different diagnostic tests among a
#' population of individuals which is partitioned into groups. The data
#' should be in the form of the output produced by \code{\link{as_epiPOMSdata}}.
#'
#' @param
#' x An object of class \sQuote{epiPOMSdata}, produced from the
#' \code{\link{as_epiPOMSdata}} function.
#' @param
#' options_gen A list of general graphical options:
#' \describe{
#' \item{main}{A string vector of length \eqn{P}, the total number of groups in
#' the study, specifying the main title for each plot. Defaults to
#' \code{sprintf("Group \%g", 1:length(} \code{x$epidatatest1))}.}
#' \item{mfrow}{A vector of two positive integers of the form \code{c(nr, nc)},
#' similar to that in \code{par}, corresponding to the number of rows and columns
#' in the resulting plotting matrix, i.e. subsequent figures will be drawn in an
#' \code{nr}-by-\code{nc} array on the device. The default
#' values are 2 for creating a 2\eqn{\times}{x}2 plotting matrix.}
#' \item{ask}{A logical value similar to that in \code{par}. If
#' \code{TRUE} (and the \R session is interactive) the user will be prompted
#' before a new page of plots is started. The default is \code{ask =}
#' \code{dev.interactive()}.}
#' }
#' @param
#' options_x A list of graphical options for the horizontal axis:
#' \describe{
#' \item{xlim}{A list containing \eqn{P} numeric vectors of length 2, giving
#' the x coordinates ranges for each plot (similar to that in
#' \code{plot}). Defaults to \code{lapply(x$epidatatest1,}
#' \code{function(x)} \code{c(1,dim(x)[2]))}, i.e. from 1 to the time of the
#' last observation in each group.}
#' \item{xlab}{A label for x-axis. Defaults to "Time".}
#' \item{xticks}{A logical value; if \code{TRUE}, tick marks on x-axis are
#' displayed at every observation time point per group. If it is \code{FALSE},
#' then the default style of the generic \code{plot} function
#' for x-axis interval calculation is used, where tick marks are equally spaced.
#' The default is \code{TRUE}.}
#' }
#' @param
#' options_y A list of graphical options for the vertical axis:
#' \describe{
#' \item{ylim}{A list containing \eqn{P} numeric vectors of length 2, giving
#' the y coordinates ranges for each plot (similar to that in
#' \code{plot}). Defaults to \code{lapply(x$epidatatest1,}
#' \code{function(x)} \code{c(1,dim(x)[1]))}, i.e. from 1 to the total number of
#' individuals per group.}
#' \item{ylab}{A label for y-axis. Defaults to "Individual".}
#' \item{yticks}{A logical value; if \code{TRUE}, tick marks on y-axis are
#' displayed at every individual index per group. If it is \code{FALSE}, then
#' the default style of the generic \code{plot} function for
#' y-axis interval calculation is used. The default is \code{TRUE}.}
#' }
#' @param
#' options_tests A list of graphical options related to test results:
#' \describe{
#' \item{cex}{A positive numerical value, similar to that in \code{plot}, giving
#'  the amount by which the observed test result labels should be magnified.
#'  Defaults to 0.6.}
#' \item{col}{A vector of two colours used for plotting the individual's observed
#' results of test 1 and test 2, respectively. By default uses red for test 1
#' and blue for test 2.}
#' }
#' @param
#' ... Additional arguments that are passed to the generic
#' \code{plot} function.
#'
#' @return
#' A series of plots, containing the observed test results for each
#' individual over the sampling time with separate plots for each group.
#' The first line for each individual represents the outcome of the first test
#' and the second line the outcome of the second test.
#'
#' @details
#' Only works for objects of class \sQuote{epiPOMSdata}. The default plotting
#' parameter values work well for epidemics up to about 10 individuals per
#' group over the sampling period of about 100 days. For a larger number of
#' individuals per group or sampling period, it is recommended to
#' change the plotting matrix values or use \code{pdf}
#' and adjust the plotting dimensions.
#'
#' @author
#' Panayiota Touloupou
#'
#' @seealso
#' \code{\link{as_epiPOMSdata}} for producing an object of class
#' \sQuote{epiPOMSdata}.
#'
#' @examples
#' # Load E. coli O157:H7 data
#' data(Ecoliepidata)
#' # Load your own dataframe using myData <- read.csv('myDataFile.csv').
#'
#' # Format the data
#' epidata <- as_epiPOMSdata(Ecoliepidata)
#'
#' # Plot the data
#' plot(epidata)
#'
#' @export
plot.epiPOMSdata <- function(x, options_gen = list(main = sprintf("Group %g", 1:length(x$epidatatest1)), 
    mfrow = c(2, 2), ask = dev.interactive()), options_x = list(xlim = lapply(x$epidatatest1, function(x) c(1, 
    dim(x)[2])), xlab = "Time", xticks = TRUE), options_y = list(ylim = lapply(x$epidatatest1, 
    function(x) c(1, dim(x)[1])), ylab = "Individual", yticks = TRUE), options_tests = list(cex = 0.6, 
    col = c("red", "blue")), ...) {
    
    # Error checks for input arguments
    if (!is(x, "epiPOMSdata")) {
        stop("plot.epiPOMSdata: The x must be a class of 'epiPOMSdata'.", call. = FALSE)
    }
    
    check <- options_gen$main
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The main at options_gen has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("plot.epiPOMSdata: The main at options_gen must be a vector.", call. = FALSE)
    } else if (length(check) != length(x$epidatatest1)) {
        stop("plot.epiPOMSdata: Length of main at options_gen is not compatible.", call. = FALSE)
    } else if (!is.character(check)) {
        stop("plot.epiPOMSdata: The entries of main at options_gen must be  characters.", call. = FALSE)
    }
    
    check <- options_gen$mfrow
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The mfrow at options_gen has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("plot.epiPOMSdata: The mfrow at options_gen must be a vector.", call. = FALSE)
    } else if (length(check) != 2) {
        stop("plot.epiPOMSdata: Length of mfrow at options_gen is not compatible.", call. = FALSE)
    } else if (suppressWarnings(any(is.na(as.numeric(check))))) {
        stop("plot.epiPOMSdata: The entries of mfrow at options_gen must be positive integers.", 
            call. = FALSE)
    } else if (any(check <= 0) | any(check != round(check))) {
        stop("plot.epiPOMSdata: The entries of mfrow at options_gen must be positive integers.", 
            call. = FALSE)
    }
    
    check <- options_gen$ask
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The ask at options_gen has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSdata: The ask at options_gen must be logical, either TRUE or FALSE.", 
            call. = FALSE)
    } else if (!is.logical(check) | is.na(check)) {
        stop("plot.epiPOMSdata: The ask at options_gen must be logical, either TRUE or FALSE.", 
            call. = FALSE)
    }
    
    check <- options_x$xlim
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The xlim at options_x has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("plot.epiPOMSdata: The xlim at options_x must be a list.", call. = FALSE)
    } else if (length(check) != length(x$epidatatest1)) {
        stop("plot.epiPOMSdata: Length of xlim at options_x is not compatible", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) ((!is.vector(x)) | is.list(x)))))) {
        stop("plot.epiPOMSdata: The xlim at options_x must contain P vectors of length 2.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (length(x) != 2))))) {
        stop("plot.epiPOMSdata: The xlim at options_x must contain P vectors of length 2.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) suppressWarnings(any(is.na(as.numeric(x)))))))) {
        stop("plot.epiPOMSdata: The entries of xlim at options_x must be real  numbers.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (sort(x) != x))))) {
        stop("plot.epiPOMSdata: The xlim at options_x must contain P vectors in ascending order.", 
            call. = FALSE)
    }
    
    check <- options_x$xlab
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The xlab at options_x has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSdata: The xlab at options_x must be a character.", call. = FALSE)
    } else if (!is.character(check)) {
        stop("plot.epiPOMSdata: The xlab at options_x must be a character.", call. = FALSE)
    }
    
    check <- options_x$xticks
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The xticks at options_x has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSdata: The xticks at options_x must be logical, either TRUE or FALSE.", 
            call. = FALSE)
    } else if (!is.logical(check) | is.na(check)) {
        stop("plot.epiPOMSdata: The xticks at options_x must be logical, either TRUE or FALSE.", 
            call. = FALSE)
    }
    
    check <- options_y$ylim
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The ylim at options_y has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("plot.epiPOMSdata: The ylim at options_y must be a list.", call. = FALSE)
    } else if (length(check) != length(x$epidatatest1)) {
        stop("plot.epiPOMSdata: Length of ylim at options_y is not compatible", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) ((!is.vector(x)) | is.list(x)))))) {
        stop("plot.epiPOMSdata: The ylim at options_y must contain P vectors of length 2.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (length(x) != 2))))) {
        stop("plot.epiPOMSdata: The ylim at options_y must contain P vectors of length 2.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) suppressWarnings(any(is.na(as.numeric(x)))))))) {
        stop("plot.epiPOMSdata: The entries of ylim at options_y must be real  numbers.", call. = FALSE)
    } else if (any(unlist(lapply(check, FUN = function(x) (sort(x) != x))))) {
        stop("plot.epiPOMSdata: The ylim at options_y must contain P vectors in ascending order.", 
            call. = FALSE)
    }
    
    check <- options_y$ylab
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The ylab at options_y has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSdata: The ylab at options_y must be a character.", call. = FALSE)
    } else if (!is.character(check)) {
        stop("plot.epiPOMSdata: The ylab at options_y must be a character.", call. = FALSE)
    }
    
    check <- options_y$yticks
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The yticks at options_y has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSdata: The yticks at options_y must be logical, either TRUE or FALSE.", 
            call. = FALSE)
    } else if (!is.logical(check) | is.na(check)) {
        stop("plot.epiPOMSdata: The yticks at options_y must be logical, either TRUE or FALSE.", 
            call. = FALSE)
    }
    
    check <- options_tests$cex
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The cex at options_tests has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSdata: The cex at options_tests must be a positive real number.", call. = FALSE)
    } else if (!is.numeric(check)) {
        stop("plot.epiPOMSdata: The cex at options_tests must be a positive real number.", call. = FALSE)
    } else if ((check <= 0)) {
        stop("plot.epiPOMSdata: The cex at options_tests must be a positive real number.", call. = FALSE)
    }
    
    
    check <- options_tests$col
    if (is.null(check)) {
        stop("plot.epiPOMSdata: The col at options_tests has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("plot.epiPOMSdata: The col at options_tests must be a vector.", call. = FALSE)
    } else if (length(check) != 2) {
        stop("plot.epiPOMSdata: Length of col at options_tests is not compatible.", call. = FALSE)
    } else if (any(!areColors(check))) {
        stop("plot.epiPOMSdata: The entries of col at options_tests must be colours.", call. = FALSE)
    }
    
    # End of error checks for input arguments
    
    
    DataTest1 <- x$epidatatest1
    DataTest2 <- x$epidatatest2
    
    ngroups <- length(DataTest1)
    
    parold <- par(no.readonly = TRUE)
    on.exit(par(parold), add = TRUE)
    par(mfrow = options_gen$mfrow)
    op <- par(ask = options_gen$ask)
    
    for (i in 1:ngroups) {
        ninds <- dim(DataTest1[[i]])[1]
        
        obsertimes <- which(!is.na(DataTest1[[i]][1, ]))
        for (j in 2:ninds) {
            obsertimes <- c(obsertimes, which(!is.na(DataTest1[[i]][j, ])))
        }
        obsertimes <- sort(unique(obsertimes))
        
        
        if ((options_y$yticks == TRUE) & (options_x$xticks == TRUE)) {
            plot(0, 0, xlim = options_x$xlim[[i]], ylim = c(options_y$ylim[[i]][1] - 0.3, options_y$ylim[[i]][2] + 
                0.3), main = options_gen$main[i], xlab = options_x$xlab, ylab = options_y$ylab, 
                xaxt = "n", yaxt = "n", ...)
            axis(side = 1, labels = obsertimes, at = obsertimes)
            axis(side = 2, labels = 1:ninds, at = 1:ninds)
        } else if ((options_y$yticks == FALSE) & (options_x$xticks == TRUE)) {
            plot(0, 0, xlim = options_x$xlim[[i]], ylim = c(options_y$ylim[[i]][1] - 0.3, options_y$ylim[[i]][2] + 
                0.3), main = options_gen$main[i], xlab = options_x$xlab, ylab = options_y$ylab, 
                xaxt = "n", ...)
            axis(side = 1, labels = obsertimes, at = obsertimes)
        } else if ((options_y$yticks == TRUE) & (options_x$xticks == FALSE)) {
            plot(0, 0, xlim = options_x$xlim[[i]], ylim = c(options_y$ylim[[i]][1] - 0.3, options_y$ylim[[i]][2] + 
                0.3), main = options_gen$main[i], xlab = options_x$xlab, ylab = options_y$ylab, 
                yaxt = "n", ...)
            axis(side = 2, labels = 1:ninds, at = 1:ninds)
        } else {
            plot(0, 0, xlim = options_x$xlim[[i]], ylim = c(options_y$ylim[[i]][1] - 0.3, options_y$ylim[[i]][2] + 
                0.3), main = options_gen$main[i], xlab = options_x$xlab, ylab = options_y$ylab, 
                ...)
        }
        
        for (j in 1:ninds) {
            Test1 <- DataTest1[[i]][j, ]
            Test2 <- DataTest2[[i]][j, ]
            pos1 <- which(!is.na(Test1))
            pos2 <- which(!is.na(Test2))
            if (length(pos1) > 0) {
                text(pos1, cex = options_tests$cex, rep(j + 0.17, length(pos1)), Test1[pos1], col = options_tests$col[1])
            }
            if (length(pos2) > 0) {
                text(pos2, cex = options_tests$cex, rep(j - 0.17, length(pos2)), Test2[pos2], col = options_tests$col[2])
            }
        }
    }
    par(op)
}

areColors <- function(x) {
    sapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
    })
}


