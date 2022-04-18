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
/* ECO387E Problem Set 3, 2 */
/* Paul Le Tran, plt377 */ 
/* 17 April, 2022 */
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
/** Loading in data and setting up variables **/
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
drop if age < 18 | age > 65
drop if class == 1 | class == 2 | class == 3 | class == 6 | class == 7

/* Returning to home directory */
cd `home_dir'
/******************************************************************************/

/******************************************************************************/
/** Part 2a: Plot the kernel densities of the predicted wage for both the **/
/** employed and unemployed on the same graph. How do you interpret the **/
/** differences in the densities? **/
/* Creating regressand */
gen ln_wages = log(wages)
/* Creating experience squared variable */
gen potexpsq = potexp^2

/* Performing WLS */
regress ln_wage schoolyrs potexp potexpsq black hispanic white ///
  i.female##i.married childpresence state* ///
  [pweight = earnwt]
  
/* Obtaining part 1b regression predictions of ln(wage) for both employed and */
/* unemployed */
predict ln_wages_hat

/* Creating variable of predicted ln(wages) for employed */
gen ln_wages_hat_emp = ln_wages_hat if lfsr94 == 1 | lfsr94 == 2

/* Creating variable of predicted ln(wages) for unemployed */
gen ln_wages_hat_unemp = ln_wages_hat if lfsr94 == 3 | lfsr94 == 4

/* Converting predictions to wages in levels */
gen wages_hat = exp(ln_wages_hat)

/* Creating variable of predicted wages for employed */
gen wages_hat_emp = wages_hat if lfsr94 == 1 | lfsr94 == 2

/* Creating variable of predicted wages for unemployed */
gen wages_hat_unemp = wages_hat if lfsr94 == 3 | lfsr94 == 4

/* Creating kernel density plots for both sets of predicted ln(wages) */
kdensity ln_wages_hat_emp, addplot(kdensity ln_wages_hat_unemp) ///
  title("Estimated kernel densities of predicted ln(wages)") ///
  xtitle("Predicted ln(wages)") ///
  legend(ring(0) pos(2) label(1 "Employed") label(2 "Unemployed"))
graph export "\\path\to\graphics\2ai_plot.png", replace
graph close

/* Creating kernel density plots for both sets of predicted wages */
kdensity wages_hat_emp, addplot(kdensity wages_hat_unemp) ///
  xtitle("Predicted Wages") ///
  title("Estimated kernel densities of predicted wages") ///
  legend(ring(0) pos(2) label(1 "Employed") label(2 "Unemployed"))
graph export "\\path\to\graphics\2aii_plot.png", replace
graph close

/**********/
/* ANSWER */
/**********/
/* In the kernel density plot of predicted ln(wages), we see that the density */
/* for those employed is a little shifted to the left. Furthermore, said */
/* density is more in-line with the shape of a normal distribution. In */
/* contrast, the density estimate for those employed has a peak more spread */
/* out and a fatter right tail. When observing the plot of predicted wages in */
/* levels, we see similar differences, though both plots are much more skewed */
/* to the right. */

/* This difference between the densities in levels and in natural logarithmic */
/* terms might have to do with the outliers affecting the distribution of the */
/* predicted wages much more so than when in log terms. As for the difference */
/* observed between groups in either plot, it seems this indicates that the */
/* predicted wages for unemployed people are much more volatile and not as */
/* "uniform" as those seen in the distribution for employed people. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 2b: Compute the average predicted wage for each type of unemployed **/
/** (job losers, temporary job ended, job leavers, re-entrants, and new **/
/** entrants). Are there any meaningful differences in predicted wages among **/
/** these groups? If yes, how do you interpret these differences? **/
/* Computing mean of predicted ln(wages) according to different unemployment */
/* types */
sum ln_wages_hat_unemp if untype == 1 | untype == 2
sum ln_wages_hat_unemp if untype == 3
sum ln_wages_hat_unemp if untype == 4
sum ln_wages_hat_unemp if untype == 5
sum ln_wages_hat_unemp if untype == 6

/**********/
/* ANSWER */
/**********/
/* When examining the predicted means of ln(wages) for each unemployment */
/* group, we see that the biggest difference is between job losers and new */
/* entrants. However, said difference is only by about 0.3. Otherwise, there */
/* really isn't any meaningful difference. */
/**************/
/* END ANSWER */
/**************/

/* Computing mean of predicted wages according to different unemployment */
/* types */
sum wages_hat_unemp if untype == 1 | untype == 2
sum wages_hat_unemp if untype == 3
sum wages_hat_unemp if untype == 4
sum wages_hat_unemp if untype == 5
sum wages_hat_unemp if untype == 6

/**********/
/* ANSWER */
/**********/
/* When examining the predicted means of wages in levels for each group, we */
/* see some meaningful differences. First off, we see a difference of almost */
/* 120 USD between job losers and new entrants, with the latter having the */
/* higher mean. This could possibly be due to new entrants being "fresh" in */
/* the labour force, so they could be viewed as having more updated skills */
/* and knowledge compared to job-losers. Another difference seen is between */
/* job losers and job leavers, with the latter having the greater mean. It's */
/* possible job leavers have a greater predicted mean due to them having more */
/* favourable opportunities in line when they quit their previous job. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/

/******************************************************************************/
/** Part 2c: Compute the average predicted wage by duration of unemployment **/
/** (0-13 weeks, 14-26 weeks, and 27+ weeks). Are there any meaningful **/
/** differences in predicted wages among these groups? If yes, how do you **/
/** interpret these differences or the lack thereof? **/
/* Computing mean of predicted ln(wages) according to different unemployment */
/* duration groups */
sum ln_wages_hat_unemp if prunedur <= 13
sum ln_wages_hat_unemp if prunedur >= 14 & prunedur <= 26
sum ln_wages_hat_unemp if prunedur >= 27

/**********/
/* ANSWER */
/**********/
/* When observing the predicted mean of ln(wages) according to different  */
/* unemployment duration by the specified groups, we see very little (if any) */
/* meaningful differences between the three averages. This could be due to */
/* natural logarithm making any variation harder to see immediately. */
/**************/
/* END ANSWER */
/**************/

/* Computing mean of predicted wages according to different unemployment */
/* duration groups */
sum wages_hat_unemp if prunedur <= 13
sum wages_hat_unemp if prunedur >= 14 & prunedur <= 26
sum wages_hat_unemp if prunedur >= 27

/**********/
/* ANSWER */
/**********/
/* When observing the predicted mean of wages according to different  */
/* unemployment duration by the specified groups, the differences become both */
/* more noticeable and interesting. Specifically, we see that the longer one */
/* is unemployed, the smaller the predicted average wage in levels. A */
/* possible explanation for this is because those whom are newly or recently */
/* unemployed still has relevant skills desired at large in the labour */
/* market. Those who are unemployed longer might be a bit more "outdated" or */
/* have other characteristics that make it difficult for them to find a job */
/* in general. */
/**************/
/* END ANSWER */
/**************/
/******************************************************************************/
