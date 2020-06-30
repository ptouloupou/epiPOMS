#' Plot the outputs of an \sQuote{epiPOMSmcmc} object.
#'
#' @description
#' The plot method of \sQuote{epiPOMSmcmc} objects (i.e. the output of the
#' \code{\link{epiPOMS_mcmc}} function) can be used to visualise two types of
#' information. The first one shows the trace and density plots for each of
#' the model parameters. The second one shows the posterior probability of
#' colonisation over time for each individual in the study.
#'
#' @param
#' x An object of class \sQuote{epiPOMSmcmc}, produced from the
#' \code{\link{epiPOMS_mcmc}} function.
#' @param
#' infoepidata A list containing general information of the epidemic data of
#' class \sQuote{epiPOMSdata}. It should be the same list used as an argument
#' of function \code{\link{epiPOMS_mcmc}}.
#' @param
#' plottype A character string specifying which of the two outputs of the
#' \code{\link{epiPOMS_mcmc}} function to plot, namely "trace" for producing
#' trace and density plots for each of the model parameters, "prob" for
#' plotting the posterior probability of colonisation over time for each
#' individual within a group or "both" for plotting both outputs. The default
#' is \code{plottype = "trace"}.
#' @param
#' colstates A vector of \eqn{(n_g + 1)} colours used for plotting the
#' probability of an individual to be at carriage state
#' \eqn{s = 0, 1, \ldots, n_g}, if \code{plottype = "prob"} or \code{"both"},
#' where \eqn{n_g} denotes the number of different strains in the study.
#' By default uses \code{colstates =}
#' \code{c("grey", rainbow(infoepidata$ntypes))}.
#' @param
#' auto.layout A logical value similar to that in \code{plot.mcmc}. If
#' \code{TRUE} automatically creates a plot layout. The default is
#' \code{TRUE}.
#' @param
#' ask A logical value similar to that in \code{par}. If
#' \code{TRUE} (and the \R session is interactive) the user will be prompted
#' before a new page of plots is started. The default is \code{ask =}
#' \code{dev.interactive()}.
#' @param
#' ... Additional arguments that are passed to the generic
#' \code{plot} function.
#'
#' @return
#' The argument \code{plottype} is used to specify what plot to produce.
#' If set to \code{"trace"}, trace and density plots of the posterior
#' distribution of the model parameters are produced using the \code{plot} S3
#' method from the \code{coda} package. If set to \code{"prob"}, the posterior
#' probabilities of an individual being colonised by a specific strain or not
#' being colonised over the sampling period are provided, with separate plots
#' for each individual within a group.  Setting \code{plottype = "both"}, both
#' outputs are plotted.
#'
#' @details
#' Using the augmented states of carriage in each recorded MCMC iteration,
#' one can estimate the probability that an individual is colonised by a
#' specific strain or not colonised for every day in the study, regardless of
#' whether diagnostic tests were collected on that day. Setting
#' \code{plottype = "prob"} or \code{"both"}, the posterior
#' probability of colonisation for each individual in the study over the
#' sampling period is provided. The plot for each individual is divided into
#' two panels. The bottom panel contains the posterior probabilities of
#' colonisation. In the top panel the observed test results are included, which
#' are imperfect measures of the true underlying process. The first line
#' represents the outcome of the samples from the first test and the second line
#' represents the outcome of second test; "-" indicates negative sample, "+"
#' indicates that the sample was positive but not chosen for strain typing,
#' otherwise, strain number is given.
#'
#' @author
#' Panayiota Touloupou
#'
#' @seealso
#' \code{\link{epiPOMS_mcmc}} for generating posterior samples of the
#' model parameters and \code{\link{summary.epiPOMSmcmc}} for
#' displaying summary information about an \sQuote{epiPOMSmcmc} object.
#'
#' @export
plot.epiPOMSmcmc <- function(x, infoepidata, plottype = "trace", colstates = c("grey", rainbow(infoepidata$ntypes)),
auto.layout = TRUE, ask = dev.interactive(), ...) {
    
    # Error checks for input arguments
    if (!is(x, "epiPOMSmcmc")) {
        stop("plot.epiPOMSmcmc: The x must be a class of 'epiPOMSmcmc'.", call. = FALSE)
    }
    
    check <- infoepidata
    checknames <- c("ntypes", "obserdata", "ngroups", "ninds", "tmax", "indNA", "tmaxNA")
    if (is.null(check)) {
        stop("plot.epiPOMSmcmc: The infoepidata has to be specified.", call. = FALSE)
    } else if (!is.list(check)) {
        stop("plot.epiPOMSmcmc: The infoepidata must be a list obtained by function info_epiPOMSdata.",
        call. = FALSE)
    } else if (length(check) != 7) {
        stop("plot.epiPOMSmcmc: The infoepidata must be a list obtained by function info_epiPOMSdata.",
        call. = FALSE)
    } else if (any(names(check) != checknames)) {
        stop("plot.epiPOMSmcmc: The infoepidata must be a list obtained by function info_epiPOMSdata.",
        call. = FALSE)
    }
    
    check <- plottype
    if (is.null(check)) {
        stop("plot.epiPOMSmcmc: The plottype has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSmcmc: The plottype must be a character string, 'trace' or 'prob' or 'both'.",
        call. = FALSE)
    } else if (!(check %in% c("trace", "prob", "both"))) {
        stop("plot.epiPOMSmcmc: The plottype must be a character string, 'trace' or 'prob' or 'both'.",
        call. = FALSE)
    }
    
    check <- colstates
    if (is.null(check)) {
        stop("plot.epiPOMSmcmc: The colstates has to be specified.", call. = FALSE)
    } else if ((!is.vector(check)) | is.list(check)) {
        stop("plot.epiPOMSmcmc: The colstates must be a vector.", call. = FALSE)
    } else if (length(colstates) != (infoepidata$ntypes + 1)) {
        stop("plot.epiPOMSmcmc: Length of colstates is not compatible.", call. = FALSE)
    } else if (any(!areColors(check))) {
        stop("plot.epiPOMSmcmc: The entries of colstates must be colours.", call. = FALSE)
    }
    
    
    check <- auto.layout
    if (is.null(check)) {
        stop("plot.epiPOMSmcmc: The auto.layout has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSmcmc: The auto.layout must be logical, either TRUE or FALSE.", call. = FALSE)
    } else if (!is.logical(check) | is.na(check)) {
        stop("plot.epiPOMSmcmc: The auto.layout must be logical, either TRUE or FALSE.", call. = FALSE)
    }
    
    check <- ask
    if (is.null(check)) {
        stop("plot.epiPOMSmcmc: The ask has to be specified.", call. = FALSE)
    } else if ((length(check) != 1) | is.list(check)) {
        stop("plot.epiPOMSmcmc: The ask must be logical, either TRUE or FALSE.", call. = FALSE)
    } else if (!is.logical(check) | is.na(check)) {
        stop("plot.epiPOMSmcmc: The ask must be logical, either TRUE or FALSE.", call. = FALSE)
    }
    
    # End of error checks for input arguments
    
    parold <- par(no.readonly = TRUE)
    on.exit(par(parold), add = TRUE)
    
    if (plottype %in% c("trace", "both")) {
        
        if (is.null(x$gammamcmc)) {
            pars <- cbind(x$alphamcmc, x$betamcmc, x$mumcmc, x$deltamcmc, x$numcmc, x$thetamcmc)
            colnames(pars) <- c(sprintf("alpha_%g", 1:infoepidata$ntypes), sprintf("beta_%g", 1:infoepidata$ntypes),
            sprintf("mu_%g", 1:infoepidata$ntypes), "delta", sprintf("nu_%g", 0:infoepidata$ntypes),
            sprintf("theta_%g", 1:2), "theta_C", "theta_S", "theta_P")
        } else {
            pars <- cbind(x$alphamcmc, x$betamcmc, x$mumcmc, x$deltamcmc, x$gammamcmc, x$numcmc,
            x$thetamcmc)
            colnames(pars) <- c(sprintf("alpha_%g", 1:infoepidata$ntypes), sprintf("beta_%g", 1:infoepidata$ntypes),
            sprintf("mu_%g", 1:infoepidata$ntypes), "delta", "gamma", sprintf("nu_%g", 0:infoepidata$ntypes),
            sprintf("theta_%g", 1:2), "theta_C", "theta_S", "theta_P")
        }
        plot(as.mcmc(pars), ...)
    }
    
    if (plottype %in% c("prob", "both")) {
        parold <- par(no.readonly = TRUE)
        
        for (pp in 1:infoepidata$ngroups) {
            
            if (auto.layout) {
                mfrow <- setmfrow(infoepidata$ninds[pp])
            }
            
            par(oma = c(2, 1, 1, 1))
            par(mfrow = mfrow)
            par(ask = ask)
            
            for (cc in 1:infoepidata$ninds[pp]) {
                counts <- x$probstates[[pp]][, cc, ]
                counts <- counts[order(nrow(counts):1), ]
                rownames(counts) <- infoepidata$ntypes:0
                colnames(counts) <- 1:infoepidata$tmax[pp]
                cols <- rev(colstates)
                xx <- barplot(counts, main = "", xlab = "Time", ylab = "Probability of colonisation",
                col = cols, legend = FALSE, space = 0, axis.lty = 1, border = NA, xlim = c(-5,
                infoepidata$tmax[pp]), ...)
                
                colstest <- c(rev(cols), "black")
                labelss <- c("-", 1:infoepidata$ntypes, "+")
                mtext(sprintf("Group %g - Individual %g", pp, cc), side = 3, line = 2.5, cex = 0.8)
                mtext("T1", side = 3, line = 1.2, at = -5, font = 2, cex = 0.6)
                mtext("T2", side = 3, line = 0.5, at = -5, font = 2, cex = 0.6)
                
                Test1 <- infoepidata$obserdata$test1[[pp]][cc, ]
                Test2 <- infoepidata$obserdata$test2[[pp]][cc, ]
                
                pos1 <- which(!is.na(Test1) & Test1 > 0 & Test1 <= infoepidata$ntypes)
                pos2 <- which(!is.na(Test2) & Test2 > 0 & Test2 <= infoepidata$ntypes)
                
                if (length(pos1) > 0) {
                    mtext(labelss[as.numeric(Test1[pos1]) + 1], col = colstest[as.numeric(Test1[pos1]) +
                    1], side = 3, line = 1.2, at = xx[pos1], font = 2, cex = 0.6)
                }
                if (length(pos2) > 0) {
                    mtext(labelss[as.numeric(Test2[pos2]) + 1], col = colstest[as.numeric(Test2[pos2]) +
                    1], side = 3, line = 0.5, at = xx[pos2], font = 2, cex = 0.6)
                }
                
                pos1 <- which(!is.na(Test1) & Test1 == 0)
                pos2 <- which(!is.na(Test2) & Test2 == 0)
                if (length(pos1) > 0) {
                    mtext(labelss[as.numeric(Test1[pos1]) + 1], col = colstest[as.numeric(Test1[pos1]) +
                    1], side = 3, line = 1.2, at = xx[pos1], font = 2, cex = 0.6)
                }
                if (length(pos2) > 0) {
                    mtext(labelss[as.numeric(Test2[pos2]) + 1], col = colstest[as.numeric(Test2[pos2]) +
                    1], side = 3, line = 0.5, at = xx[pos2], font = 2, cex = 0.6)
                }
                
                pos1 <- which(!is.na(Test1) & Test1 == (infoepidata$ntypes + 1))
                pos2 <- which(!is.na(Test2) & Test2 == (infoepidata$ntypes + 1))
                if (length(pos1) > 0) {
                    mtext(labelss[as.numeric(Test1[pos1]) + 1], col = colstest[as.numeric(Test1[pos1]) +
                    1], side = 3, line = 1.2, at = xx[pos1], font = 2, cex = 0.6)
                }
                if (length(pos2) > 0) {
                    mtext(labelss[as.numeric(Test2[pos2]) + 1], col = colstest[as.numeric(Test2[pos2]) +
                    1], side = 3, line = 0.5, at = xx[pos2], font = 2, cex = 0.6)
                }
            }
            par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
            legend("bottom", legend = rownames(counts), fill = cols, bty = "n", horiz = TRUE, xpd = TRUE,
            inset = c(0, 0))
            
            par(parold)
        }
    }
}

areColors <- function(x) {
    sapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
    })
}

setmfrow <- function(Nplots) {
    
    if (Nplots == 1) {
        mfrow <- c(1, 1)
    } else if (Nplots == 2) {
        mfrow <- c(2, 1)
    } else if (any(Nplots == c(3, 4))) {
        mfrow <- c(2, 2)
    } else if (any(Nplots == c(5, 6, 10, 11, 12))) {
        mfrow <- c(3, 2)
    } else {
        mfrow <- c(3, 3)
    }
    return(mfrow)
}



