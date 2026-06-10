# Microwave-free vector magnetometry and crystal orientation determination with NV centers

This repository contains the MATLAB scripts used in the manuscript

**"Microwave-free vector magnetometry and crystal orientation determination with Nitrogen-Vacancy centers using Bayesian inference"**

by H. Espinós, O. Dhungel, A. Wickenbrock, D. Budker, R. Puebla, and E. Torrontegui.

The code implements the Bayesian inference framework used to reconstruct crystal orientation and magnetic-field parameters from photoluminescence (PL) maps obtained from microwave-free NV-center cross-relaxation measurements.

## Overview

The repository contains research scripts developed specifically for the analyses reported in the manuscript. They are provided to facilitate reproducibility of the published results and are not intended as a general-purpose software package.

The implementation consists of:

* A forward model describing photoluminescence changes induced by NV–NV cross-relaxation resonances.
* A Gaussian likelihood function including uncertainty propagation from experimental control variables.
* A brute-force Bayesian parameter estimation procedure based on grid evaluation of the posterior distribution.
* Computation of marginal posterior distributions for the inferred parameters.

## Files

### `log_likelihood_derivatives.m`

Evaluates the log-likelihood of experimental photoluminescence data given a set of model parameters.

Main features:

* Analytical model of cross-relaxation resonances.
* Lorentzian photoluminescence features.
* Numerical evaluation of derivatives with respect to experimental control variables.
* Propagation of uncertainties in magnetic-field and angular coordinates.
* Gaussian likelihood evaluation.

### `brute_force_general.m`

Performs Bayesian parameter estimation through brute-force evaluation of the likelihood on a multidimensional parameter grid.

Main tasks:

* Loading experimental photoluminescence data.
* Defining parameter grids.
* Parallel evaluation of the likelihood function.
* Computation of posterior distributions.
* Marginalization over nuisance parameters.
* Visualization of one- and two-dimensional posterior distributions.

## Requirements

The code was developed and tested using:

* MATLAB R2025a
* Parallel Computing Toolbox (recommended for efficient execution)

## Data

Experimental data used in the publication are available through the Zenodo repository associated with the manuscript.

The scripts expect photoluminescence maps in CSV format. File paths should be adapted to the user's local directory structure.

## Reproducibility

The scripts reproduce the Bayesian analysis framework presented in the manuscript. Some experimental parameters (e.g., calibration factors, noise estimates, parameter ranges, and dataset paths) are specified directly in the scripts and should be adjusted according to the dataset being analyzed.
