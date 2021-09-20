

---------------------------------------------------------------------------------------------------------------------------
This is the README.txt file for R scripts to implement analyses in paper:
'Estimating the Optimal Timing of Surgery from Observational Data'.


Date: 2020-01-18
---------------------------------------------------------------------------------------------------------------------------

R scripts:

1. DesignMatrixFun.R: the function for constructing spline function design matrix A.
2. IPSWEst.R: calculating the generalized propensity score with competing-risk model for main analysis.
3. IPSWEstSensitivity.R: calculating the generalized propensity score with competing-risk model for sensitivity analysis.
4. WeightedEst.R: Cox ph model with IPSW and model selection. 

How to run:

1. Get the public SVRT data set from NIH/NHLBI Pediatric Heart Network website:
   https://www.pediatricheartnetwork.com/pud_login.asp?study_id=SVR.

2. Use the IPSWEst.R to create the generalized propensity score.

3. Fit the Cox ph model with IPSW by WeightedEst.R.

4. For the sensitivity analysis, which has competing risks elective S2P, non-elective S2P, death/heart transplantation,
   please run IPSWEstSensitivity.R to get the generalized propensity score. 

---------------------------------------------------------------------------------------------------------------------------

