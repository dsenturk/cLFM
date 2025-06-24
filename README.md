# cLVM

## Contents of this folder

- cLFM_Tutorial.Rmd: A step-by-step implementation of cLFM and associated procedures including data generating process of the first simulation scenario and evaluation of model parameter estimation. Details of simulation design, EM estimation of model parameters and trimming-refinement process are described in Section 3 and Section 4 of 'Contrastive Latent Functional Model'.
- Simulation_dgp.R: Functions for generating contrastive functional data sets with latent components and score variances of Scenario 1 described in the simulation section in 'Contrastive Latent Functional Model'.
- cLFM_Functions_1.R: Functions for conducting the proposed EM estimation algorithm described in Algorithm 1 of 'Contrastive Latent Functional Model'.
- cLFM_Functions_2.R: Functions for implementing cLFM modeling with the trimming-refinement process described in Algorithm 2 of 'Contrastive Latent Functional Model'.

## Introduction

The contents of this folder allow for the implementation of cLFM for contrastive functional data settings as proposed in "Contrastive Latent Functional Model". Users can simulate contrastive data pairs as described in Simulation Scenario 1 and apply the proposed estimation algorithm to fit the cLFM. Detailed instructions on how to perform the aforementioned procedures are included in cLFM_Tutorial.Rmd.

## Requirements

- The included R programs require R 4.4.2 and the packages and files listed in cLFM_Tutorial.Rmd.

## Installation

Load the R program files into the global environment and install the required packages using commands in cLFM_Tutorial.Rmd.
