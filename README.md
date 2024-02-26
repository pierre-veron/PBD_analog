# PBD_analog
Test the analogy between the protracted birth-death model and the birth-death 
model

## Overview
Simulations and statistical inferences of phylogenies under the protracted 
birth-death (PBD) model following the works from 
> Etienne, R. S., & Rosindell, J. (2012). Prolonging the Past Counteracts the Pull of the Present: Protracted Speciation Can Explain Observed Slowdowns in Diversification. _Systematic Biology, 61_, 204.

> Etienne, R. S., Morlon, H., & Lambert, A. (2014). Estimating the duration of speciation from phylogenies. _Evolution, 68_, 2430–2440.

Test of the analog expressions of birth and death rates from the PBD model. 

## Project architecture
* `PBD_bayesian_fit_BD` : scripts to fit a birth-death model to trees simulated from the PBD model, with Bayesian approach.
* `fig` : figures of the project
* `PBD_notebook.ipynb` : worksheet for the tests of the equivalent diverfication rates.
* `modules/PBD_analog.py` : implementation of the equivalent birth-death rates from the PBD model and function to simulate a PBD process. 
* `test_predictions`: test empirically the predictions derived from PBD model (without trees so far, only simple stochastic processes). 

## Contact
* Pierre Veron (`pierre.veron.2017@polytechnique.org`)
* Jérémy Andréoletti