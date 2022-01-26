# Analysis of Critical Values of Response Rate (RR) in Single-Arm Two-Stage Phase II Clinical Trials

Yi Liu 

November 2021

In this project, we investigated the pattern of testing boundary of respond rate (RR) in two-stage phase II trials. Although the boundary can be controlled by some relationships
that have been known in large sample theory, we also want to see their trends under small samples, because phase II trials require small samples mainly for ethical consideration.
The results of this project to some degree can help us to find some new criteria for making choice between minimax and optimal design of single-arm two-stage phase II trials.

Slides:
- RR_Analysis_Ph2Trials.pdf: Contains the statistical background, method introduction, scientific question of interest, numerical experiments and results, and conclusion. 

Code and data: 
- design_FS_func.R: Function for single-arm two-stage phase II trials design (minimax and optimal) with only futility (lower) stopping
- design_twoS_func.R: Function for single-arm two-stage phase II trials design (minimax and optimal) with both futility (lower) and superiority (upper) stopping
- data_FS.R: Generate data for designs (minimax and optimal) with only futility stopping by different choices of parameters (can be found in slides)
- data_twoS.R: Generate data for designs (minimax and optimal) with both futility and superiority stopping by different choices of parameters (can be found in slides)
- design_FS.Rdata: Generated data by data_FS.R, the parameters used are those in slides
- design_twoS.Rdata: Generated data by data_twoS.R, the parameters used are those in slides
- plot_FS.R: Generate trend plots and summary statistics by data from data_FS.R
- plot_twoS.R: Generate trend plots and summary statistics by data from data_twoS.R
