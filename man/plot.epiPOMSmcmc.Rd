% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.epiPOMSmcmc.R
\name{plot.epiPOMSmcmc}
\alias{plot.epiPOMSmcmc}
\title{Plot the outputs of an \sQuote{epiPOMSmcmc} object.}
\usage{
\method{plot}{epiPOMSmcmc}(x, infoepidata, plottype = "trace",
  colstates = c("grey", rainbow(infoepidata$ntypes)), auto.layout = TRUE,
  ask = dev.interactive(), ...)
}
\arguments{
\item{x}{An object of class \sQuote{epiPOMSmcmc}, produced from the
\code{\link{epiPOMS_mcmc}} function.}

\item{infoepidata}{A list containing general information of the epidemic data of
class \sQuote{epiPOMSdata}. It should be the same list used as an argument
of function \code{\link{epiPOMS_mcmc}}.}

\item{plottype}{A character string specifying which of the two outputs of the
\code{\link{epiPOMS_mcmc}} function to plot, namely "trace" for producing
trace and density plots for each of the model parameters, "prob" for
plotting the posterior probability of colonisation over time for each
individual within a group or "both" for plotting both outputs. The default
is \code{plottype = "trace"}.}

\item{colstates}{A vector of \eqn{(n_g + 1)} colours used for plotting the
probability of an individual to be at carriage state
\eqn{s = 0, 1, \ldots, n_g}, if \code{plottype = "prob"} or \code{"both"},
where \eqn{n_g} denotes the number of different strains in the study.
By default uses \code{colstates =}
\code{c("grey", rainbow(infoepidata$ntypes))}.}

\item{auto.layout}{A logical value similar to that in \code{plot.mcmc}. If
\code{TRUE} automatically creates a plot layout. The default is
\code{TRUE}.}

\item{ask}{A logical value similar to that in \code{par}. If
\code{TRUE} (and the \R session is interactive) the user will be prompted
before a new page of plots is started. The default is \code{ask =}
\code{dev.interactive()}.}

\item{...}{Additional arguments that are passed to the generic
\code{plot} function.}
}
\value{
The argument \code{plottype} is used to specify what plot to produce.
If set to \code{"trace"}, trace and density plots of the posterior
distribution of the model parameters are produced using the \code{plot} S3
method from the \code{coda} package. If set to \code{"prob"}, the posterior
probabilities of an individual being colonised by a specific strain or not
being colonised over the sampling period are provided, with separate plots
for each individual within a group.  Setting \code{plottype = "both"}, both
outputs are plotted.
}
\description{
The plot method of \sQuote{epiPOMSmcmc} objects (i.e. the output of the
\code{\link{epiPOMS_mcmc}} function) can be used to visualise two types of
information. The first one shows the trace and density plots for each of
the model parameters. The second one shows the posterior probability of
colonisation over time for each individual in the study.
}
\details{
Using the augmented states of carriage in each recorded MCMC iteration,
one can estimate the probability that an individual is colonised by a
specific strain or not colonised for every day in the study, regardless of
whether diagnostic tests were collected on that day. Setting
\code{plottype = "prob"} or \code{"both"}, the posterior
probability of colonisation for each individual in the study over the
sampling period is provided. The plot for each individual is divided into
two panels. The bottom panel contains the posterior probabilities of
colonisation. In the top panel the observed test results are included, which
are imperfect measures of the true underlying process. The first line
represents the outcome of the samples from the first test and the second line
represents the outcome of second test; "-" indicates negative sample, "+"
indicates that the sample was positive but not chosen for strain typing,
otherwise, strain number is given.
}
\seealso{
\code{\link{epiPOMS_mcmc}} for generating posterior samples of the
model parameters and \code{\link{summary.epiPOMSmcmc}} for
displaying summary information about an \sQuote{epiPOMSmcmc} object.
}
\author{
Panayiota Touloupou
}
