We use the following files in obtaining the simulation results reported in the main paper.  All results were obtained in Matlab.

Data:

FRQMainMat.csv 	- Balkrishnan et al. (2014) data

To reproduce the simulation results, one would run files in the following order after appropriately modifying the files to refer to directories on the user's computer:

1.  FRQFitModel.m                  - main .m file for fitting spatio-temporal model to Balakrishnan et al. (2014) data
2.  FRQFitMean.m                   - get OLS estimates of regression coefficients, fixed effects, and parameters on functions of lagged errors
3.  FRQSim_GenY.m                  - file that generates simulated Y values
4.  FRQSim_EstAR.m                 - estimate ad hoc "autoregressive" models in simulated data for use in ad hoc adaptive procedure
5.  FRQSim_Main.m                  - obtain most simulation results
6.  FRQSim_FinalStats.m            - calculate most summaries of simulation results reported in paper and supplement, includes all calculation of Canay-Romano-Shaikh with two year time block groups
7.  FRQSim_WildBoot.m              - obtain wild bootstrap critical values  
8.  FRQSim_WildBootStats.m         - calculate size using wild bootstrap critical values
9.  FRQSim_CanayRomanoShaikh.m     - more Canay-Romano-Shaikh results
10. FRQSim_FinalTWFirmXYear.m      - simulation to get two-way-clustered by firm and year results
11. FRQSim_FinalStatsTWFirmXYear.m - calculate summaries of results with two-way-clustered by firm and year standard errors
12. FRQSim_MBB_RMSE1.m             - calculate RMSE with firm and time fixed effects using balanced sample used in block bootstrap results 
13. FRQSim_MBB_RMSE2.m             - calculate RMSE with firm and state x time fixed effects, firm and one digit SIC x year fixed effects, firm and two digit SIC x year fixed effects, and firm and size category x year fixed fixed effects using balanced sample used in block bootstrap results 
14. FRQSim_MBB_RMSE3.m             - calculate RMSE with firm x 8 year and year fixed effects, firm x 6 year and year fixed effects, firm x 4 year and year fixed effects, and firm x 2 year and year fixed effects using balanced sample used in block bootstrap results 
15. FRQSim_MBB1.m                  - obtain block bootstrap critical value for rows 1-5 in Table 3
16. FRQSim_MBB2.m                  - obtain block bootstrap critical value for rows 6-9 in Table 3
17. FRQSim_MBB3.m                  - obtain block bootstrap critical value for rows 10-13 in Table 3
18. FRQSim_MBB4.m                  - obtain block bootstrap critical value for rows 14-17 in Table 3
19. FRQSim_MBB1_WrapUp.m           - tabulate results using block bootstrap critical value for rows 1-5 in Table 3
20. FRQSim_MBB2_WrapUp.m           - tabulate results using block bootstrap critical value for rows 6-9 in Table 3
21. FRQSim_MBB3_WrapUp.m           - tabulate results using block bootstrap critical value for rows 10-13 in Table 3
22. FRQSim_MBB4_WrapUp.m           - tabulate results using block bootstrap critical value for rows 14-17 in Table 3
23. FRQSim_MBB_Power.m             - calculate power using MBB critical values.  
24. FRQSim_FinalStats_Gauss.m      - calculate size using standard normal critical values, results reported in supplementary appendix


Dependent data inference files:
cluster_se.m		- calculate clustered standard errors
FamaMacbeth.m		- perform Fama-MacBeth estimation
FamaMacbethRand.m	- Canay-Romano-Shaikh permutation inference
cluster_wildboot.m	- wild bootstrap with Mammen weights critical values
index_mbb.m		- generate indices for blocks to use in panel overlapping block bootstrap


Auxiliary Files:
FRQSim_GenErr.m		- called by FRQSim_GenY.m to generate error terms for constructing simulated outcomes
FRQFillLagGroupMean.m	- generate functions of lagged errors used in defining lagged part of covariance structure
FRQFillLagGroupMeanTT.m - generate functions of lagged errors used in defining lagged part of covariance structure
FRQLagGroupMean.m	- generate functions of lagged errors used in defining lagged part of covariance structure
groupcross.m		- create variable where entries correspond to categories formed by fully interacting two discrete variables
recode.m		- helper file that recodes discrete variable with K unique entries to run from 1:K
egb2_pdf.m		- pdf of EGB2 distribution
egb2rnd.m		- generate draws from EGB2 random variable
NumToLatex.m            - convert number to latex friendly form


Matlab Data Files used in generating simulated outcomes:
bcv_optim1.mat-bcv_optim6.mat 	- intermediate files with optimization results produced by FRQFitModel.m
distparams.mat			- estimated parameters for EGB2 distribution produced by FRQFitModel.m
SigmaEsts.mat			- estimated covariance matrices produced FRQFitModel.m
XCatIndices.mat			- indices of categories used to form lagged variables produced by FRQFitMean.m
TimeIndices.mat			- time indices produced by FRQFitMean.m
lagparams.mat			- estimated parameters for "autoregressive" part of model produced by FRQFitMean.m
fixedeffects.mat		- estimates firm and time fixed effects produced by FRQFitMean.m
FirmIndices.mat			- Firm indices produced by FRQFitMean.m
xMat.mat			- Design matrix (X) produced by FRQFitMean.m
olsparams.mat			- OLS estimates of parameters on X produced by FRQFitMean.m


Matlab Data Files generated in running the simulation:
BCVResidualRegressions.mat			- Estimated parameters used in ad hoc adaptive procedure produced by FRQSim_EstAR.m
FRQSim_FINAL_PART1.mat				- Most simulation results.  Produced by FRQSim_Main.m
FRQSim_FINAL_PART2.mat				- Most simulation results including tabulation of reported summary statistics
FRQSim_FINAL_BOOT.mat				- Bootstrap critical values
FRQSim_FINAL_FIRMTIME.mat			- Simulation results using two-way clustered by firm and time standard errors
FRQSim_FINAL_FIRMTIME1.mat      		- Simulation results including tabulation of two-way clustered by firm and time standard errors
FRQSim_FINAL_PART2_GAUSS.mat			- Most simulation results including tabulation of reported summary statistics using Gaussian critical values
SizePt1_1000X2000_bsize3_17Obs_sim_1000.mat 	- block bootstrap critical value for rows 1-5 in Table 3
SizePt2_1000X2000_bsize3_17Obs_sim_1000.mat 	- block bootstrap critical value for rows 6-9 in Table 3
SizePt3_1000X2000_bsize3_17Obs_sim_1000.mat 	- block bootstrap critical value for rows 10-13 in Table 3
SizePt4_1000X2000_bsize3_17Obs_sim_1000.mat 	- block bootstrap critical value for rows 14-17 in Table 3
