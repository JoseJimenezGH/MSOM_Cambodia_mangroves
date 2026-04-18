R code implementing a Bayesian multi‑species occupancy model (Site × Year × Session × Species) for camera‑trap data from 2017–2023 using JAGS.
Occupancy and detection are modelled hierarchically with species and year random effects, species‑specific temporal trends, environmental covariates, and habitat effects.
Data are stacked only where sampling effort is positive, ensuring inference is conditioned on observed survey effort.
