################
##### INIT #####
################

#Load packages
require(car) #linearHypothesis
require(estimatr) #lm_robust and iv_robust
require(haven) #read_dta
require(marginaleffects) #avg_comparisons
require(MatchIt) #matchit
require(optmatch) #matchit
require(sandwich) #avg_comparisons
require(tidyverse) #tidy

#Setting random seed for replicability
set.seed(12345)

#Set working directory
#setwd("your path here")

#Read in the data
BCG14 = read_dta("replication_data.dta")

#################################
##### VARIABLE DESCRIPTIONS #####
#################################

# CLA: Indicator of the person being assigned to the standard track (control group)
# CVE: Indicator of the person being assigned to the public new track
# OPP: Indicator of the person being assigned to the private new track
# EMPLOI_6MOIS: Indicator of exiting unemployment to employment after 6 months
# POIDS_PZ_6MOIS: Inverse of the product of estimated assignement probabilities
# acceptationCVE_6MOIS: Indicator of the person following the public new track
# acceptationOPP_6MOIS: Indicator of the person following the private new track

### Some covariates
# woman: Indicator of the person being female
# nivetude2: Indicator of the person having finished highschool
# marie: Indicator of the person being married
# nochild: Indicator of the person having no children
# PersLayoff: Indicator of the person having experienced a personal layoff before entering unemployment
# PersLayoff: Indicator of the person being French

#Vector of all covariates
X = c("nivetude1", "nivetude3", "nivetude4", "Cadre", "Techn", "EmployQ", "EmployNQ", "OuvrQ", 
       "agegr2635", "agegr3645", "agegr4655", "agegr56", "femme", "marie", "onechild", "twoormorechild", 
       "French", "African", "IdF", "North", "ce1", "ce2", "EconLayoff", "PersLayoff", "EndCDD", 
       "EndInterim", "exper0", "exper1_5", "rsqstat2", "rsqstat3", "tempcomp", "dezus", "salaireB", 
       "salaireC", "salaireD", "salaireE", "salaireG", "primo", "Insertion", "Interim", "Q1", "Q2", "Q3")

#Vector of selected covariates
Xall = c("nivetude2", "woman", "marie", "nochild", "PersLayoff", "French")

#####################
##### BALANCING #####
#####################
for (i in Xall) {
  balancing = lm_robust(as.formula(paste(i, "~ 0 + CLA + CVE + OPP")), 
                            se_type = "stata", data = BCG14)
  #Display estimates
  print(tidy(balancing, conf.int = TRUE))
  #Test if levels of CLA and CVE are equal
  print(linearHypothesis(balancing, "CLA = CVE"))
  #Test if levels of CLA and OPP are equal
  print(linearHypothesis(balancing, "CLA = OPP"))
  #Test if levels of CVE and OPP are equal
  print(linearHypothesis(balancing, "CVE = OPP"))
}

######################
##### COMPLIANCE #####
######################

#Zero terms suppress constants

#Public program (Table 4, column 1, scaled with factor 100) 
#Estimate the model
compliance_m1 = lm_robust(acceptationCVE_6MOIS ~ 0 + CVE + OPP + CLA, weights = POIDS_PZ_6MOIS, 
                           se_type = "stata", data = BCG14)
#Display estimates
tidy(compliance_m1, conf.int = TRUE)

#Private program (Table 4, column 3, scaled with factor 100) 
#Estimate the model
compliance_m2 = lm_robust(acceptationOPP_6MOIS ~ 0 + CVE + OPP + CLA, weights = POIDS_PZ_6MOIS, 
                           se_type = "stata", data = BCG14)
#Display estimates
tidy(compliance_m2, conf.int = TRUE)

######################################
##### INTENTION-TO-TREAT EFFECTS #####
######################################

#Non-weighted OLS regression
#Estimate the model
ITT_m1 = lm_robust(EMPLOI_6MOIS ~ CVE + OPP, se_type = "stata", data = BCG14)
#Display estimates
tidy(ITT_m1, conf.int = TRUE)

#Weighted OLS regression (Table 3, column 1, scaled with factor 100) 
#Estimate the model
ITT_m2 = lm_robust(EMPLOI_6MOIS ~ CVE + OPP, weights = POIDS_PZ_6MOIS, se_type = "stata", data = BCG14)
#Display estimates
tidy(ITT_m2, conf.int = TRUE)
#Test if ITT estimates of CVE and OPP are equal
linearHypothesis(ITT_m2, "CVE = OPP")

#Weighted OLS regression with controls (Table 3, column 2, scaled with factor 100) 
#Write formula
formula_ITT_m3 = as.formula(paste("EMPLOI_6MOIS ~ CVE + OPP + ", paste(X, collapse= "+")))
#Display formula
formula_ITT_m3
#Estimate the model
ITT_m3 = lm_robust(formula_ITT_m3, weights = POIDS_PZ_6MOIS, se_type = "stata", data = BCG14)
#Display estimates
tidy(ITT_m3, conf.int = TRUE)
#Test if ITT estimates of CVE and OPP are equal
linearHypothesis(ITT_m3, "CVE = OPP")

###########################################
##### LOCAL AVERAGE TREATMENT EFFECTS #####
###########################################

#General note: IV SEs in R are slightly larger than those in Stata

#Weighted IV regression (Table 5, column 1, scaled with factor 100) 
#Estimate the model
LATE_m1 = iv_robust(EMPLOI_6MOIS ~ acceptationCVE_6MOIS + acceptationOPP_6MOIS | CVE + OPP,
                     weights = POIDS_PZ_6MOIS, se_type = "stata", data = BCG14)
#Display estimates
tidy(LATE_m1, conf.int = TRUE)
#Test if ITT estimates of CVE and OPP are equal
linearHypothesis(LATE_m1, "acceptationCVE_6MOIS = acceptationOPP_6MOIS")

#Weighted IV regression with controls (Table 5, column 2, scaled with factor 100) 
#Write formula 
formula_LATE_m2 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationCVE_6MOIS + acceptationOPP_6MOIS + ", 
                       paste(X, collapse= "+")), paste("CVE + OPP + ", paste(X, collapse= "+")), 
                       sep = " | "))
#Display formula
formula_LATE_m2
#Estimate the model
LATE_m2 = iv_robust(formula_LATE_m2, weights = POIDS_PZ_6MOIS, se_type = "stata", data = BCG14)
#Display estimates
tidy(LATE_m2, conf.int = TRUE)
#Test if ITT estimates of CVE and OPP are equal
linearHypothesis(LATE_m2, "acceptationCVE_6MOIS = acceptationOPP_6MOIS")

#################################################
##### CONDITIONAL AVERAGE TREATMENT EFFECTS #####
#################################################

#LATE for women (Table 6, column 1, scaled with factor 100)
#Write formula 
formula_CATE_m1 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationCVE_6MOIS + acceptationOPP_6MOIS + ", 
                                   paste(X, collapse= "+")), paste("CVE + OPP + ", paste(X, collapse= "+")), 
                                   sep = " | "))
#Display formula
formula_CATE_m1
#Estimate the model
CATE_m1 = iv_robust(formula_CATE_m1, weights = POIDS_PZ_6MOIS, se_type = "stata", 
                    data = BCG14[which(BCG14$woman == 1), ])
#Display estimates
tidy(CATE_m1, conf.int = TRUE)
#Test if ITT estimates of CVE and OPP are equal
linearHypothesis(CATE_m1, "acceptationCVE_6MOIS = acceptationOPP_6MOIS")

#LATE for men (Table 6, column 2, scaled with factor 100)
#Write formula 
formula_CATE_m2 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationCVE_6MOIS + acceptationOPP_6MOIS + ", 
                                   paste(X, collapse= "+")), paste("CVE + OPP + ", paste(X, collapse= "+")), 
                                   sep = " | "))
#Display formula
formula_CATE_m2
#Estimate the model
CATE_m2 = iv_robust(formula_CATE_m1, weights = POIDS_PZ_6MOIS, se_type = "stata", 
                    data = BCG14[which(BCG14$woman == 0), ])
#Display estimates
tidy(CATE_m2, conf.int = TRUE)
#Test if ITT estimates of CVE and OPP are equal
linearHypothesis(CATE_m2, "acceptationCVE_6MOIS = acceptationOPP_6MOIS", singular.ok = TRUE)

#####################################
##### PROPENSITY SCORE MATCHING #####
#####################################

#Keeping a random sample of 10,000 observations to save time
BCG14 <- BCG14[sample(1:43977, 10000), ]

#ATE of public track on new employment
#Write matching formula
formula_PSM_matchit1 = as.formula(paste("acceptationCVE_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_PSM_matchit1
#Estimate propensity scores
PSM_matchit1 <- matchit(formula_PSM_matchit1, data = BCG14[which(BCG14$acceptationOPP_6MOIS != 1), ],
              method = "full", estimand = "ATE")
#Write equation formula
formula_PSM_m1 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationCVE_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_PSM_m1
#Estimate the model
PSM_m1 = lm(formula_PSM_m1, weights = weights, data = match.data(PSM_matchit1))
#Estimate ATE
avg_comparisons(PSM_m1, variables = "acceptationCVE_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(PSM_matchit1), acceptationCVE_6MOIS == 1), wts = "weights")

#ATT/ATET of public track on new employment
#Write matching formula
formula_PSM_matchit2 = as.formula(paste("acceptationCVE_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_PSM_matchit2
#Estimate propensity scores
PSM_matchit2 <- matchit(formula_PSM_matchit2, data = BCG14[which(BCG14$acceptationOPP_6MOIS != 1), ],
                        method = "full", estimand = "ATT")
#Write equation formula
formula_PSM_m2 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationCVE_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_PSM_m2
#Estimate the model
PSM_m2 = lm(formula_PSM_m2, weights = weights, data = match.data(PSM_matchit2))
#Estimate ATE
avg_comparisons(PSM_m2, variables = "acceptationCVE_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(PSM_matchit2), acceptationCVE_6MOIS == 1), wts = "weights")

#ATE of private track on new employment
#Write matching formula
formula_PSM_matchit3 = as.formula(paste("acceptationOPP_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_PSM_matchit3
#Estimate propensity scores
PSM_matchit3 <- matchit(formula_PSM_matchit3, data = BCG14[which(BCG14$acceptationCVE_6MOIS != 1), ],
                        method = "full", estimand = "ATE")
#Write equation formula
formula_PSM_m3 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationOPP_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_PSM_m3
#Estimate the model
PSM_m3 = lm(formula_PSM_m3, weights = weights, data = match.data(PSM_matchit3))
#Estimate ATE
avg_comparisons(PSM_m3, variables = "acceptationOPP_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(PSM_matchit3), acceptationOPP_6MOIS == 1), wts = "weights")

#ATT/ATET of private track on new employment
#Write matching formula
formula_PSM_matchit4 = as.formula(paste("acceptationOPP_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_PSM_matchit4
#Estimate propensity scores
PSM_matchit4 <- matchit(formula_PSM_matchit4, data = BCG14[which(BCG14$acceptationCVE_6MOIS != 1), ],
                        method = "full", estimand = "ATT")
#Write equation formula
formula_PSM_m4 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationOPP_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_PSM_m4
#Estimate the model
PSM_m4 = lm(formula_PSM_m4, weights = weights, data = match.data(PSM_matchit4))
#Estimate ATE
avg_comparisons(PSM_m4, variables = "acceptationOPP_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(PSM_matchit4), acceptationOPP_6MOIS == 1), wts = "weights")

#####################################
##### NEAREST-NEIGHBOR MATCHING #####
#####################################

#ATC of public track on new employment
#Write matching formula
formula_NN_matchit1 = as.formula(paste("acceptationCVE_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_NN_matchit1
#Estimate propensity scores
NN_matchit1 <- matchit(formula_NN_matchit1, data = BCG14[which(BCG14$acceptationOPP_6MOIS != 1), ],
                        estimand = "ATC", distance = "mahalanobis")
#Write equation formula
formula_NN_m1 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationCVE_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_NN_m1
#Estimate the model
NN_m1 = lm(formula_NN_m1, weights = weights, data = match.data(NN_matchit1))
#Estimate ATC
avg_comparisons(NN_m1, variables = "acceptationCVE_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(NN_matchit1), acceptationCVE_6MOIS == 1), wts = "weights")

#ATT/ATET of public track on new employment
#Write matching formula
formula_NN_matchit2 = as.formula(paste("acceptationCVE_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_NN_matchit2
#Estimate propensity scores
NN_matchit2 <- matchit(formula_NN_matchit2, data = BCG14[which(BCG14$acceptationOPP_6MOIS != 1), ],
                        estimand = "ATT", distance = "mahalanobis")
#Write equation formula
formula_NN_m2 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationCVE_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_NN_m2
#Estimate the model
NN_m2 = lm(formula_NN_m2, weights = weights, data = match.data(NN_matchit2))
#Estimate ATE
avg_comparisons(NN_m2, variables = "acceptationCVE_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(NN_matchit2), acceptationCVE_6MOIS == 1), wts = "weights")

#ATC of private track on new employment
#Write matching formula
formula_NN_matchit3 = as.formula(paste("acceptationOPP_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_NN_matchit3
#Estimate propensity scores
NN_matchit3 <- matchit(formula_NN_matchit3, data = BCG14[which(BCG14$acceptationCVE_6MOIS != 1), ],
                        estimand = "ATC", distance = "mahalanobis")
#Write equation formula
formula_NN_m3 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationOPP_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_NN_m3
#Estimate the model
NN_m3 = lm(formula_NN_m3, weights = weights, data = match.data(NN_matchit3))
#Estimate ATC
avg_comparisons(NN_m3, variables = "acceptationOPP_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(NN_matchit3), acceptationOPP_6MOIS == 1), wts = "weights")

#ATT/ATET of private track on new employment
#Write matching formula
formula_NN_matchit4 = as.formula(paste("acceptationOPP_6MOIS ~", paste(X, collapse= "+")))
#Display formula
formula_NN_matchit4
#Estimate propensity scores
NN_matchit4 <- matchit(formula_NN_matchit4, data = BCG14[which(BCG14$acceptationCVE_6MOIS != 1), ],
                        estimand = "ATT", distance = "mahalanobis")
#Write equation formula
formula_NN_m4 = as.formula(paste(paste("EMPLOI_6MOIS ~ acceptationOPP_6MOIS*(", paste(X, collapse= "+")), ")"))
#Display formula
formula_NN_m4
#Estimate the model
NN_m4 = lm(formula_NN_m4, weights = weights, data = match.data(NN_matchit4))
#Estimate ATC
avg_comparisons(NN_m4, variables = "acceptationOPP_6MOIS", vcov = ~subclass,
                newdata = subset(match.data(NN_matchit4), acceptationOPP_6MOIS == 1), wts = "weights")

