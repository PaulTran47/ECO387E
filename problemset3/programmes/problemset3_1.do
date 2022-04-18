/******************************************************************************/
/* 1 asterisk represent commented out code */
/** 2 asterisks represent comments for code or commented out code **/

/* One asterisk represents commented out code. */
/** Two asterisks represent comments for code. **/

/* Answers to question parts that don't involve code can be found at the */
/* bottom of the programme, in the section ``Questions asked in problemset */
/* that don't involve code. */

/* Text answers to question parts that involve code will be between the */
/* sub-section label: */
/**********/
/* ANSWER */
/**********/
/* Answer here */
/**************/
/* END ANSWER */
/**************/

/* Comments that are important will be between the sub-section label: */
/********/
/* NOTE */
/********/
/* Important note here */
/************/
/* END NOTE */
/************/
/* ECO387E Problem Set 3, 1 */
/* Paul Le Tran, plt377 */ 
/* 14 April, 2022 */
/******************************************************************************/

/******************************************************************************/
/** Setting up workspace **/
clear
cls

set matsize 547

local home_dir = "\\path\to\programmes"
local data_dir = "\\path\to\data"

cd `home_dir'
/******************************************************************************/

/******************************************************************************/
/** Part 1a: Loading in data and setting up variables **/
cd `data_dir'
use morg18

/* Creating dummy variables as specified by the problem set */
/* Female dummy */
gen female = 0
replace female = 1 if sex == 2

/* Hispanic dummy (not caring about nationality */
gen hispanic = 0
replace hispanic = 1 if !missing(ethnic)

/* Black dummy */
gen black = 0
replace black = 1 if race == 2

/* White dummy */
gen white = 0
replace white = 1 if race == 1

/* Married dummy */
gen married = 0
replace married = 1 if marital == 1 | marital == 2 | marital == 3

/* Presence of at least one child dummy */
gen childpresence = 0
replace childpresence = 1 if ch02 == 1 | ch05 == 1 | ch35 == 1 | ch613 == 1 | ch1417 == 1

/* State dummy */
tab stfips, generate(state)

/* Occuptation dummy */
tab docc00, generate(occupation)

/* Variable showing years of education */
gen schoolyrs = 0
replace schoolyrs = 0 if ((grade92 == 31) - (grade92 == 00))
replace schoolyrs = 3 if (grade92 == 32)
replace schoolyrs = 7 if ((grade92 == 33) - (grade92 == 34))
replace schoolyrs = 9 if (grade92 == 35)
replace schoolyrs = 10 if (grade92 == 36)
replace schoolyrs = 11 if (grade92 == 37)
replace schoolyrs = 12 if (grade92 == 38)
replace schoolyrs = 12 if (grade92 == 39)
replace schoolyrs = 13.5 if (grade92 == 40)
replace schoolyrs = 14 if ((grade92 == 41) - (grade92 == 42))
replace schoolyrs = 16 if (grade92 == 43)
replace schoolyrs = 17.5 if (grade92 == 44)
replace schoolyrs = 17.5 if (grade92 == 45)
replace schoolyrs = 18 if (grade92 == 46)

/* Creating variable showing potential experience */
gen potexp = age - schoolyrs - 6

/* Creating and winsorising hourly wages variable */
/********/
/* NOTE */
/********/
/* For workers not paid hourly wages, we will find the imputed hourly wages */
/* by dividing weekly earnings with hours worked last week. */
/************/
/* END NOTE */
/************/
gen wages = 0
replace wages = earnhre if paidhre == 1
replace wages = earnwke/hourslw if paidhre == 2
replace wages =. if wages == 0
/* Windsorising hourly earnings */
/********/
/* NOTE */
/********/
/* At the moment, we will be setting the values below and above the 10th and */
/* 90th percentiles, respecitively, to the value specified by the problem */
/* set. */
/************/
/* END NOTE */
/************/
sum earnwke
local earnwke_max = r(max)
sum wages, detail
local wages_p10 = r(p10)
local wages_p90 = r(p90)
replace wages = (`earnwke_max'*1.4/35) if wages < `wages_p10' | wages > `wages_p90'

/* Dropping observations correlating to full-time students, those ages 17 and */
/* and younger, those ages 66 and older, government employees, self-employed */
/* works, and self-incorporated workers from the sample. */
drop if studftpt == 1
drop if age <= 17 | age >= 66
drop if class == 1 | class == 2 | class == 3 | class == 6 | class == 7

/* Returning to home directory */
cd `home_dir'
/******************************************************************************/

/******************************************************************************/
/** Part 1b: Carry out and show the results of a linear weighted regression **/
/** of the log hourly wage on years of schooling, experience, experience **/
/** squared, dummies for black, white, female, married, married*female, **/
/** presence of children, and state dummies. Be sure to use the earnings **/
/** weights as the survey weight in the regression. What is the R2 of the **/
/** regression? How do you interpret the dummy on marital status for both **/
/** men and women? **/

/* Creating regressand */
gen ln_wages = log(wages)
/* Creating experience squared variable */
gen potexpsq = potexp^2

/* Performing WLS */
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* ///
  [pweight = earnwt]

/**********/
/* ANSWER */
/**********/
/* The R^{2} of the regression is 0.1075. The coefficient associated with the */
/* marital status dummy can be interpreted to mean that for males with zero */
/* years of schooling and zero units of experience, being married is */
/* associated with a wage decrease of about 17.9219% (-0.179219 is the */
/* y-intercept of ln(wage)). However, this coefficient on the married dummy */
/* isn't the full effect on wages for females. This is because for married */
/* females with zero schooling and experience, they should see an associated */
/* decrease in wages by -17.9219 + 4.2784 + 06.2090 = 7.4345%. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 1c: Carry out and show the results of a linear weighted regression **/
/** with the same explanatory variables as in part 1b, but adding dummies **/
/** for each occupation. How does the R2 change relative to part 1b? How **/
/** does the coefficient on years of schooling change relative to part 1b? **/
/** Do you have an explanation for the latter? **/
/* Performing WLS */
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* occupation* ///
  [pweight = earnwt]

/**********/
/* ANSWER */
/**********/
/* The R^{2} of the regression increases to 0.2499 when compared to that of */
/* part 1b. Furthermore, we see that the coefficient on years of schooling */
/* increases from -0.1457 to -0.1417. It should be mentioned that this change */
/* isn't much at all in terms of magnitude. It’s possible this (albeit small) */
/* increase is due to how when accounting for both experience and occupation, */
/* the effect of schooling on one’s wage is more “powerful”. However, it’s */
/* also entirely possible that the observed change is small to begin with */
/* because omitted variables is still explaining much of the variation that */
/* the regression can’t account for, and therefore weakens any possible */
/* effect education has on one’s wages. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 1d: Compute three different versions of the Mean-min ratio for the **/
/** raw hourly wage, where you use the raw hourly wage, where you use **/
/** percentiles 1, 5, and 10 for the min. Compute the standard deviation of **/
/** the log hourly wage. **/
sum wages, detail
local wages_p1 = r(p1)
local wages_p5 = r(p5)
local wages_p10 = r(p10)
local wages_mean = r(mean)
gen wages_meanmin_p1 = `wages_mean'/`wages_p1'
gen wages_meanmin_p5 = `wages_mean'/`wages_p5'
gen wages_meanmin_p10 = `wages_mean'/`wages_p10'

sum ln_wages
local ln_wages_sd = r(sd)
/**********/
/* ANSWER */
/**********/
/* Using percentiles 1, 5, and 10 for the min of the raw hourly wages, we get */
/* a mean-min ratio of 28.5612, 22.5113, and 16.7084 respectively. The */
/* standard deviation of the log hourly wage is 1.4539. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 1e: For both regressions in parts 1b and 1c, compute the residuals **/
/** of the regression and compute the standard deviation of the residuals. **/
/** Then, add the sample mean of the log hourly wage to the residuals and **/
/** take the exponent of it. Compute three different versions of the **/
/** Mean-min ratio based on these residualised wages, where you use **/
/** percentiles 1, 5, and 10 for the min. **/
/****************************/
/** For part 1b regression **/
/****************************/
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* ///
  [pweight = earnwt]

/* Obtaining regression residuals */
predict resid_1b, residuals

/* Adding sample mean of ln(wages) to the residuals and taking the exponent */
sum ln_wages
gen resid_wages_1b = exp(resid_1b + r(mean))

/* Calculating the Mean-min ratios */
sum resid_wages_1b, detail
local resid_wages_1b_p1 = r(p1)
local resid_wages_1b_p5 = r(p5)
local resid_wages_1b_p10 = r(p10)
local resid_wages_1b_mean = r(mean)
gen resid_wages_1b_meanmin_p1 = `resid_wages_1b_mean'/`resid_wages_1b_p1'
gen resid_wages_1b_meanmin_p5 = `resid_wages_1b_mean'/`resid_wages_1b_p5'
gen resid_wages_1b_meanmin_p10 = `resid_wages_1b_mean'/`resid_wages_1b_p10'
/**********/
/* ANSWER */
/**********/
/* Using percentiles 1, 5, and 10 for the min of the residualised hourly */
/* wages, we get mean-min ratios of 30.8904, 18.2498, and 12.3754 */
/* respectively. */
/**************/
/* END ANSWER */
/**************/

/****************************/
/** For part 1c regression **/
/****************************/
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* occupation* ///
  [pweight = earnwt]
  
/* Obtaining regression residuals */
predict resid_1c, residuals

/* Adding sample mean of ln(wages) to the residuals and taking the exponent */
sum ln_wages
gen resid_wages_1c = exp(resid_1c + r(mean))

/* Calculating the Mean-min ratios */
sum resid_wages_1c, detail
local resid_wages_1c_p1 = r(p1)
local resid_wages_1c_p5 = r(p5)
local resid_wages_1c_p10 = r(p10)
local resid_wages_1c_mean = r(mean)
gen resid_wages_1c_meanmin_p1 = `resid_wages_1c_mean'/`resid_wages_1c_p1'
gen resid_wages_1c_meanmin_p5 = `resid_wages_1c_mean'/`resid_wages_1c_p5'
gen resid_wages_1c_meanmin_p10 = `resid_wages_1c_mean'/`resid_wages_1c_p10'
/**********/
/* ANSWER */
/**********/
/* Using percentiles 1, 5, and 10 for the min of the residualised hourly */
/* wages, we get mean-min ratios of 50.1957, 23.8760, and 14.9134 */
/* respectively. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 1f: Repeat the steps 1b - 1e, but exclude all observations where **/
/** the hourly wages are below the federal hourly wage of $7.25/hour. How **/
/** sensitive are the distributional statistics in parts 1d and 1e to **/
/** excluding these observations? **/
/* Creating and winsorising hourly wages variable */
/********/
/* NOTE */
/********/
/* For workers not paid hourly wages, we will find the imputed hourly wages */
/* by dividing weekly earnings with hours worked last week. */
/************/
/* END NOTE */
/************/
drop wages
gen wages = 0
replace wages = earnhre if paidhre == 1
replace wages = earnwke/hourslw if paidhre == 2
replace wages =. if wages == 0
/* Dropping all wages below the federal hourly wage */
drop if wages < 7.25
/* Windsorising hourly earnings */
/********/
/* NOTE */
/********/
/* At the moment, we will be setting the values below and above the 10th and */
/* 90th percentiles, respecitively, to the value specified by the problem */
/* set. */
/************/
/* END NOTE */
/************/
sum earnwke
local earnwke_max = r(max)
sum wages, detail
local wages_p10 = r(p10)
local wages_p90 = r(p90)
replace wages = (`earnwke_max'*1.4/35) if wages < `wages_p10' | wages > `wages_p90'

/*******************/
/* Redoing part 1b */
/*******************/
/* Creating regressand */
drop ln_wages
gen ln_wages = log(wages)
/* Creating experience squared variable */
drop potexpsq
gen potexpsq = potexp^2

/* Performing WLS */
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* ///
  [pweight = earnwt]
  
/*******************/
/* Redoing part 1c */
/*******************/
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* occupation* ///
  [pweight = earnwt]
  
/*******************/
/* Redoing part 1d */
/*******************/
sum wages, detail
local wages_p1 = r(p1)
local wages_p5 = r(p5)
local wages_p10 = r(p10)
local wages_mean = r(mean)
gen wages_meanmin_p1_1f = `wages_mean'/`wages_p1'
gen wages_meanmin_p5_1f = `wages_mean'/`wages_p5'
gen wages_meanmin_p10_1f = `wages_mean'/`wages_p10'

sum ln_wages
local ln_wages_sd_1f = r(sd)
/**********/
/* ANSWER */
/**********/
/* Using percentiles 1, 5, and 10 for the min of the raw hourly wages we get */
/* mean-min ratios of 27.0717, 21.3857, and 15.9543 respectively. Compared */
/* to the mean-min ratios calculated using all wages (including those below */
/* the federal hourly wage), we see a difference of a little more than 1 for */
/* all three ratios. When comparing the distributional statistics, we see */
/* little differences between using the full sample of wages or only those */
/* equal to or above the federal hourly wage. */
/**************/
/* END ANSWER */
/**************/

/*******************/
/* Redoing part 1e */
/*******************/
/****************************/
/** For part 1b regression **/
/****************************/
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* ///
  [pweight = earnwt]

/* Obtaining regression residuals */
predict resid_1b1f, residuals

/* Adding sample mean of ln(wages) to the residuals and taking the exponent */
sum ln_wages
gen resid_wages_1b1f = exp(resid_1b1f + r(mean))

/* Calculating the Mean-min ratios */
sum resid_wages_1b1f, detail
local resid_wages_1b1f_p1 = r(p1)
local resid_wages_1b1f_p5 = r(p5)
local resid_wages_1b1f_p10 = r(p10)
local resid_wages_1b1f_mean = r(mean)
gen resid_wages_1b1f_meanmin_p1 = `resid_wages_1b1f_mean'/`resid_wages_1b1f_p1'
gen resid_wages_1b1f_meanmin_p5 = `resid_wages_1b1f_mean'/`resid_wages_1b1f_p5'
gen resid_wages_1b1f_meanmin_p10 = `resid_wages_1b1f_mean'/`resid_wages_1b1f_p10'
/**********/
/* ANSWER */
/**********/
/* Using percentiles 1, 5, and 10 for the min of the residualised hourly */
/* wages, we get mean-min ratios of 29.0474, 17.3144, and 11.8859 */
/* respectively. As seen in the comparison with results from part 1d, the */
/* differences between the ratios here and those of part 1e are only by about */
/* one, which isn't very much. */
/**************/
/* END ANSWER */
/**************/

/****************************/
/** For part 1c regression **/
/****************************/
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* occupation* ///
  [pweight = earnwt]
  
/* Obtaining regression residuals */
predict resid_1c1f, residuals

/* Adding sample mean of ln(wages) to the residuals and taking the exponent */
sum ln_wages
gen resid_wages_1c1f = exp(resid_1c1f + r(mean))

/* Calculating the Mean-min ratios */
sum resid_wages_1c1f, detail
local resid_wages_1c1f_p1 = r(p1)
local resid_wages_1c1f_p5 = r(p5)
local resid_wages_1c1f_p10 = r(p10)
local resid_wages_1c1f_mean = r(mean)
gen resid_wages_1c1f_meanmin_p1 = `resid_wages_1c1f_mean'/`resid_wages_1c1f_p1'
gen resid_wages_1c1f_meanmin_p5 = `resid_wages_1c1f_mean'/`resid_wages_1c1f_p5'
gen resid_wages_1c1f_meanmin_p10 = `resid_wages_1c1f_mean'/`resid_wages_1c1f_p10'
/**********/
/* ANSWER */
/**********/
/* Using percentiles 1, 5, and 10 for the min of the residualised hourly */
/* wages, we get mean-min ratios of 46.8794, 22.1818, and 14.1877 */
/* respectively. For percentiles 5 and 10, the differences between the ratios */
/* here and those of part 1e are quite small (about one or less). However, we */
/* see that the difference between the ratios when using percentile 1 as the */
/* min is by almost four. This is by far the largest difference observed. */
/**************/
/* END ANSWER */
/**************/

/**********/
/* ANSWER */
/**********/
/* Overall, we notice that both the distributional statistics in both parts */
/* 1d and 1e are not very sensitive to excluding observations whose hourly */
/* wages are below the federal hourly wage. */
/**************/
/* END ANSWER */
/**************/