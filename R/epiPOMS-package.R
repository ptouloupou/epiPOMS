#' Bayesian inference for partially observed multi-strain epidemics.
#'
#' @description
#' The \R package \pkg{epiPOMS} provides tools for inference on epidemiological
#' data using partially observed multi-strain (\acronym{POMS}) epidemic models,
#' focusing  on applications where observations are gathered longitudinally and
#' the population under investigation is organised in small groups. These models
#' are also known as coupled hidden Markov models, where the coupling between
#' different chains accounts for the interaction between individuals within a
#' group.
#'
#' The package can be used for simulating from, and performing Bayesian
#' MCMC-based inference for individual-level multi-state epidemics with partial
#' observations. The model allows for both imperfect diagnostic tests and
#' strain misclassification (in the sense that the procedure used to classify
#' strains may indicate carriage by the wrong strain), as well as
#' between strain competition. An overview of the implemented model is given by
#' \insertCite{Touloupou2020;textual}{epiPOMS}. The package also provides
#' facilities for plotting and extracting information from the data.
#'
#' @name epiPOMS-package
#'
#' @aliases epiPOMS epiPOMS-package
#'
#' @docType package
#'
#' @details
#' The key functions for this package are:
#' \describe{
#' \item{\code{\link{obserdata_sim}}}{Simulates epidemics for \acronym{POMS}
#' models.}
#' \item{\code{\link{epiPOMS_mcmc}}}{Performs Bayesian inference on parameters
#' for \acronym{POMS} epidemic models.}
#' \item{\code{\link{plot.epiPOMSmcmc}}}{Displays diagnostic plots.}
#' }
#'
#' @author
#' Panayiota Touloupou \cr
#' Maintainer: \packageMaintainer{epiPOMS}
#'
#' @references
#' \insertAllCited{}
#'
#' @keywords package
#'
#' @useDynLib epiPOMS, .registration=TRUE
#' @import methods
#' @import stats
#' @importFrom utils packageDescription
#' @importFrom Rdpack reprompt
#' @importFrom graphics axis par plot text barplot legend mtext
#' @importFrom parallel detectCores
#' @importFrom grDevices col2rgb dev.interactive rainbow
#' @importFrom foreach foreach %dopar%
#' @importFrom coda as.mcmc
NULL
