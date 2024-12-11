# Quantifying-Individual-Risk-for-Binary-Outcome
Reproducible code for simulation and application in the article "Quantifying Individual Risk for Binary Outcome".


**For Application:**
- rhc.csv: this dataset contains the right heart catheterization data used in the application.
  - The dataset is available at this link: https://hbiostat.org/data/
  - Detailed descriptions of the variables are available at this link: https://hbiostat.org/data/repo/rhc. 

- Main.R: this file contains the code for data analysis.
  - To reproduce all the results presented in the application, simply run this file. 


**For Simulation:** 

- EstFun.R: functions for estimating FNA. There are two functions and we provide a detailed description of each below:
  - EstFun(): It is designed to estimate the FNA without applying variable selection during nuisance parameter estimation.
  - EstFun2(): It is designed to estimate the FNA while incorporating variable selection for nuisance parameter estimation.


- GenData.R: data-generating mechanisms for cases (C1)-(C6). 


- Simulation codes for estimating ATE in the target population, which use "EstFun.R" to estimate the estimand and "GenData.R" to generate the simulated datasets.  
  - MainCase1.R: main code for simulation cases (C1).   
  - MainCase2.R: main code for simulation cases (C2).
  - MainCase3.R: main code for simulation cases (C3).   
  - MainCase4.R: main code for simulation cases (C4).
  - MainCase5.R: main code for simulation cases (C5).
  - MainCase6.R: main code for simulation cases (C6).

- To reproduce all simulation results, run the file prefixed with "Main" one by one.
