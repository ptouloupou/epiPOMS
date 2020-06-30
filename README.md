# epiPOMS
## Bayesian inference for partially observed multi-strain epidemics

A collection of tools for simulating and analysing partially observed multi-strain (POMS) epidemic data. The epidemic models considered are discrete-time individual based multi-state models describing the transmission dynamics of an infectious disease among a population of individuals partitioned into groups. The R package provides tools for performing Bayesian inference in this model class.

#### Description

The R package epiPOMS provides tools for inference on epidemiological data using partially observed multi-strain (POMS) epidemic models, focusing on applications where observations are gathered longitudinally and the population under investigation is organised in small groups. These models are also known as coupled hidden Markov models, where the coupling between different chains accounts for the interaction between individuals within a group.

The package can be used for simulating from, and performing Bayesian MCMC-based inference for individual-level multi-state epidemics with partial observations. The model allows for both imperfect diagnostic tests and strain misclassification (in the sense that the procedure used to classify strains may indicate carriage by the wrong strain), as well as between strain competition. An overview of the implemented model is given by Touloupou et al. (2020). The package also provides facilities for plotting and extracting information from the data.

#### Details
The key functions for this package are:

- obserdata_sim: Simulates epidemics for POMS models.

- epiPOMS_mcmc: Performs Bayesian inference on parameters for POMS epidemic models.

- plot.epiPOMSmcmc: Displays diagnostic plots.

#### References
Touloupou P, Finkenstädt Rand B, Besser TE, French NP, Spencer SEF (2020). “Bayesian Inference for multi-strain epidemics with application to *Escherichia Coli O157:H7* in feedlot cattle.” The Annals of Applied Statistics (in press).

