# Vrije Universiteit Amsterdam
# Summer School Course: Statistical Methods for Causal Inference
# Week 1, Seminar 1:
# Teacher: dr. Sanchayan Banerjee
# Teaching assistant: Marjan Nikoloski

# The dataset needed for the exercises in this seminar is used by Alan et al. (2021)
# in their paper "Building Social Cohesion in Ethnically Mixed Schools: An Intervention on Perspective Taking" and has been adjusted for our needs

# Note: This do file should be run using Stata

# Arrangements: EVV (UGlobe) Nov 2021; Marjan Nikoloski July 2023***
library(readxl)
library(tidyverse)
Alan2021 = read_excel("./Alan2021.xlsx")

# Part 1: 

# Exercise 1: Visualization of the relationship between two variables in the dataset
# 1a: scatterplot
library(ggplot2)
ggplot(Alan2021, aes(x = srefshare, y = b_schoolsize)) + geom_point() + geom_smooth(method = "lm", se = FALSE)
# Exercise 2: Running a simple linear regression, uncovering the relationship between the level of empathy of the students (dependent variable), and the class size (independent variable)
simplemodel <- lm(fempathy_pt ~ f_csize, data = Alan2021)
summary(simplemodel)
# Exercise 3: Running a multivariate regression model, with 5 additionall variables
# 3a: regression and significance
multiplemodel <- lm(fempathy_pt ~ astudent + ageinm + male + refugee + b_schoolsize + f_csize, data = Alan2021)
summary(multiplemodel)
# 3c: Variance Inflation Factor (Multicollinearity)
library(car)
vif(multiplemodel)
# 3d: Breusch-Pagan/Cook-Weisberg test (Heteroskedasticity)
library(lmtest)
bptest(multiplemodel)
library(skedastic)
cook_weisberg(multiplemodel)
# Exercise 4: Running a binary logistic model
# 4a: regression and significance
logisticmodel <- glm(fsbully_c ~ ageinm + male + srefshare + b_schoolsize, data = Alan2021, family = binomial)
summary(logisticmodel)
odds_ageinm = exp(logisticmodel$coefficients[2])
odds_male = exp(logisticmodel$coefficients[3])
odds_srefshare = exp(logisticmodel$coefficients[4])
odds_schoolsize = exp(logisticmodel$coefficients[5])
# 4b: Marginal effects of continuous variable “age in months”
library(margins)
margins(logisticmodel)

# Part 2: Randomised experiments – Simulations
# sample size of 5000
N <- 5000
# continuous variable, drawn from a normal distribution with mean 1 and s.d. 2
x1 <- rnorm(N, mean = 1, sd = 2)
# coefficients determining data generating process for outcome variable
b0 <- 5
b1 <- 0.5
# outcome variable, linear function of coefficients, x1, and standard normal
# distributed error term
y0 <- b0 + b1 * x1 + rnorm(N)
# generate potential outcomes representing a positive treamtent effect that is
# highly statistically significant
y1 <- y0 + mean(y0) + rnorm(N)
#t-test between the two outcome variables
t.test(y1, y0, paired=TRUE)
# put into dataframe
library(magrittr)
sim_data <- data.frame(y0, y1, x1)
# true value of the average treatment effect
true_treatment_effect <- mean(sim_data$y1)-mean(sim_data$y0)
true_treatment_effect
# Non-randomised treatment: 
# What if we have treatment assignment (d1) where cases with lower x1 are more likely to receive treatment?
x1_median <- median(sim_data$x1)
sim_data$d1 <- 0
sim_data$d1 <- ifelse(sim_data$x1 < x1_median, 1, sim_data$d1)
sim_data$outcome <- ifelse(sim_data$d1 == 0, sim_data$y0, sim_data$y1)
# average treatment effect
treatment_effect_1 <- lm(outcome ~ d1, data = sim_data)
summary(treatment_effect_1)
# Randomization: now we randomly assign cases to treatment
sim_data$d2 <- runif(nrow(sim_data)) < 0.5
sim_data$outcome2 <- ifelse(sim_data$d2 == 0, sim_data$y0, sim_data$y1)
# average treatment effect
treatment_effect_randomised <- lm(outcome2 ~ d2, data = sim_data)
summary(treatment_effect_randomised)

# Part 3:
# Exercise 1: Baseline balance checks for differences in background characteristics
# Exogenous
ttest_result_1 <- t.test(bnode_in_friend ~ treatment, data = Alan2021)
ttest_result_1
ttest_result_2 <- t.test(bnode_in_supportself ~ treatment, data = Alan2021)
ttest_result_2
ttest_result_3 <- t.test(bnode_in_studyself ~ treatment, data = Alan2021)
ttest_result_3
# Endogenous
ttest_result_4 <- t.test(bnode_in_studyself ~ refugee, data = Alan2021)
ttest_result_4
ttest_result_5 <- t.test(bnode_in_supportself ~ refugee, data = Alan2021)
ttest_result_5
ttest_result_6 <- t.test(bnode_in_studyself ~ refugee, data = Alan2021)
ttest_result_6
ttest_result_7 <- t.test(bnode_in_studyself ~ astudent, data = Alan2021)
ttest_result_7
ttest_result_8 <- t.test(bnode_in_supportself ~ astudent, data = Alan2021)
ttest_result_8
ttest_result_9 <- t.test(bnode_in_studyself ~ astudent, data = Alan2021)
ttest_result_9
# Exercise 2: Treatment Effects on Student and Teacher Reports of Violence and Antisocial Behavior
# In-class bullying
in_class_bully <- lm(fsbully_c ~ treatment + bsbully_c + ageinm + male + refugee + astudent + braven_sd + beyes_sd + b_schoolsize + f_csize + dist1 + dist2 + dist3 + dist4 + dist5 + dist6 + dist7 + dist8 + dist9 + dist10 + factor(bstrata), data = Alan2021, cluster = Alan2021$b_schoolid)
summary(in_class_bully)
# Out-Class Bullying
out_class_bully <- lm(fsbully_s ~ treatment + bsbully_c + ageinm + male + refugee + astudent + braven_sd + beyes_sd + b_schoolsize + f_csize + dist1 + dist2 + dist3 + dist4 + dist5 + dist6 + dist7 + dist8 + dist9 + dist10 + factor(bstrata), data = Alan2021, cluster = Alan2021$b_schoolid)
summary(out_class_bully)
# Teacher Behavioral Grade
teacher_bully <- lm(tbully_f ~ treatment + bsbully_c + ageinm + male + refugee + astudent + braven_sd + beyes_sd + b_schoolsize + f_csize + dist1 + dist2 + dist3 + dist4 + dist5 + dist6 + dist7 + dist8 + dist9 + dist10 + factor(bstrata), data = Alan2021, cluster = Alan2021$b_schoolid)
summary(teacher_bully)

# Exercise 3: Heterogeneous Treatment Effects on Student and Teacher Reports of Violence and Antisocial Behavior
# In-class bullying
in_class_bully_het <- lm(fsbully_c ~ bsbully_c + treatment + refugee + treatment:refugee, data = Alan2021)
summary(in_class_bully_het)
# Out-Class Bullying
out_class_bully_het <- lm(fsbully_s ~ bsbully_c + treatment + refugee + treatment:refugee, data = Alan2021)
summary(out_class_bully_het)
# Teacher Behavioral Grade
teacher_bully_het <- lm(tbully_f ~ bsbully_c + treatment + refugee + treatment:refugee, data = Alan2021)
summary(teacher_bully_het)