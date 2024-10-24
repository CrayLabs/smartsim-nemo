# SmartSim-NEMO

## Purpose 
This repository is a work-in-progress to recreate the MOM6 ML-EKE 
(Partee et al. 2022) work in NEMO 5.0. This takes advantage of the
recent inclusion of the GEOMETRIC framework. Much like the ML-EKE
paper, the prognostic equation is replaced by the a prediction
from a CNN using features calculated from resolved model quantities.

## Status
Initial implementation in progress. Features are being calculated
according to NEMO conventions and in a similar way to the MOM6
variables. Additionally, a new module (dbclient.F90) handles
the initialization of the SmartRedis client (called from a modified
nemogcm.F90).

## TODO
- Sanity check the calculated features
- Create a driver script for NEMO's GYRE_PISCES at 1/4-degree
- Demonstrate inference-in-the-loop and compare results to GEOMETRIC
