# SCR_trend_models
Multisession closed population Bayesian spatial capture-recapture for estimating density over time.

This repo contains an R script and source code for fitting models in parallel in Nimble, including customised configuration of the Nimble MCMC. The data formatting script is specifically designed for Panthera's Leopard Monitoring Program datasets, however those are not currently publicly available.

The source code includes three specifications of the density process: 
    1) A static density model that assumes density is constant across all sessions and that any differences in observed patterns are the result of binomial process variance or sampling variance.
    2) A random walk ("rw") density model that assumes density at time t varies randomly from density at time t-1. This is essentially equivalent to an independent density process across sessions but allows for projections into the future.
    3) A trend model that assumes density exhibits a log-linear pattern over time.

All models use the same detection process:
- Detection is modeled as a binomial process with daily occasions and detection probability decaying according to a half-normal function. This specification leverages efficient vectorized distributions from the *nimbleSCR* package.
- Sex is treated as a partially observed covariate.
- Sigma is modeled independently between the sexes.
- Baseline detection probability is modeled using a random trap effect and a centered sex effect.
- The statespace is continuous but bounded. It is square and the entire state space is assumed to be suitable habitat.

These models build on the trend models presented in Chapter 4 of Rogan ([2021](https://open.uct.ac.za/items/9ca2ca34-e0e9-4ec4-99f7-9f3151da5215)): The application of spatial capture-recapture models to investigate leopard ecology and
conservation in South Africa, a Doctoral Thesis from the University of Cape Town.

Development of this repo is ongoing.

