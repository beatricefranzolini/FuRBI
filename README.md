# FuRBI
Code used for the paper "Nonparametric priors with full-range borrowing of information" by F. Ascolani, B. Franzolini, A. Lijoi, and I. Prünster
The documentation of the various files reads:

- "Simulation_FirstPart.R" contains the code to replicate Figures 1, S2, S3 and obtain Table 1.
- "Functions.R" contains the auxiliary functions for "Simulation_FirstPart.R": in particular it contains the code for the various methods used in the Illustration.
- "Simulation_SecondPart.R" contains the code to replicate Figure 2.
- "FunctionsMVN_rho.R" contains the auxiliary functions for "Simulation_SecondPart.R".
- "Simulation_correlation.R" contains the code to replicate Figure S1.
- "LSBP.R" contains the code to replicate Figure S4

- "Furbi_GM_forMultiNA.R" contains the auxiliary functions used for estimating the model in Section 6.4 
- "main_multiNAsimul.R" contains the code to replicate the simulation studies in Section 6.4 
- "main_brandsma.R" contains the code to replicate the clustering analysis in Section 6.4 for the brandsma dataset

- "acrmgm.py" contains python functions to estimate Additive n-FuRBI mixture model.
- "data.xlsx" contains stocks and commodities data.
- "CPO.py" contains the code to compute CPO measures.
- "empirical_corr.py" contains the code to get empirical correlation plot. 
- "main.py" contains the code to estimate the additive n-FuRBI mixture model on financial data.
- "mainExch_Ind.py" contains the code to estimate exchangeable and independent mixture models on financial data.
