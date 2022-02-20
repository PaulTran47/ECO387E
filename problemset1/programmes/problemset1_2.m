%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% Answers to question parts that don't involve code can be found at the
% bottom of the programme, in the section ``Questions asked in problemset x
% that don't involve code".

% Text answers to question parts that involve code will be between the
% sub-section label:
%=======
% ANSWER
%=======
% Answer here
%===========
% END ANSWER
%===========

% Comments that are important will be between the sub-section label:
%=====
% NOTE
%=====
% Important note here
%=========
% END NOTE
%=========
% ECO387E Problem Set 1, 2
% Paul Le Tran, plt377
% 16 February, 2022
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';

cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 2a: Model with exogenous separations. Compute the steady-state (SS) values of the wage, w, and job-finding rate, \phi.
% Creating calibration variables
p_a = 1;
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding the SS unemployment rate given theta_ss (we are given it is 0.05)
u_ss_a = 0.05;

% Solving for theta that yields unemployment rate of 0.05
theta0 = 0.10;
theta_for_u005_ss = fsolve(@(theta) u_theta_ss(theta) - u_ss_a, theta0);

% Creating SS function k(\theta). This is done by substituting the wage
% equation into the job creation condition equation, then rearranging to
% have k as a function of \theta.
k_theta_ss = @(theta) ((1 - beta)*(p_a - b).*q_theta(theta))./(lambda + r + beta.*phi_theta(theta));

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p_a + theta.*k_theta_ss(theta));

% Finding SS k given the \theta value that results in SS unemployment rate
% to equal 0.05
k_ss_a = k_theta_ss(theta_for_u005_ss);

% Finding SS wage given the \theta value that results in SS unemployment
% rate to equal 0.05
w_ss_a = w_theta_ss(theta_for_u005_ss);

% Finding SS job-finding rate given the \theta value that results in SS
% unemployment rate to equal 0.05
phi_ss_a = phi_theta(theta_for_u005_ss);
%=======
% ANSWER
%=======
% Given our calibration values and mismatching parameter \mu = 1, we find
% that k = 1.6043 will result in SS equilibrium u = 0.05, \theta = 0.1343, w = 0.9871, \phi(\theta) = 0.57.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 2b: Compute the SS unemployment rate, job-finding rate, and wage for the same calibration as in 2a., but choosing p = 1.1. What is the SS elasticity of u, w, and \phi with respect to productivity p?
% Creating calibration variables
p = 1.1;
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = 0.47;
theta_ss = fsolve(@(theta) jc_theta_ss(theta), theta0);

% Finding the SS wage given theta_ss
w_ss = w_theta_ss(theta_ss);

% Finding the SS unemployment rate given theta_ss
u_ss = u_theta_ss(theta_ss);

% Finding the SS job-finding rate given theta_ss
phi_ss = phi_theta(theta_ss);
%=======
% ANSWER
%=======
% With productivity p = 1.1 and every other calibration parameter set the
% same as in part 2a., we find that SS u = 0.0479, \phi(\theta) = 0.5957, w = 1.0856.
%===========
% END ANSWER
%===========

% To avoid more mindless algebra and derivatives (ASK PROFESSOR MUELLER IF
% THIS IS ACTUALLY THE DESIRED METHOD), we will use the ratio of per cent
% changes definition of elasticity with respect to p.
u_elasticity_wrt_p = ((u_ss - u_ss_a)/(u_ss + u_ss_a))/((p - p_a)/(p + p_a));
w_elasticity_wrt_p = ((w_ss - w_ss_a)/(w_ss + w_ss_a))/((p - p_a)/(p + p_a));
phi_elasticity_wrt_p = ((phi_ss - phi_ss_a)/(phi_ss + phi_ss_a))/((p - p_a)/(p + p_a));
%=======
% ANSWER
%=======
% Using p \in {1, 1.1}, the SS elasticities of unemployment rate, wage,
% and job-finding rate with respect to p are -0.4402, 0.9975, and 0.4628
% respectively.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 2c: Compute the SS elasticities of u, w, and \phi with respect to p for the same calibration as in 2a., but where you set \beta = 0.2. Moreover, is the SS unemployment rate with this alternative calibration too low or too high to the socially efficient unemployment rate? Describe your intuition why this is so.
%======================================================
% Calculating SS u, w, and \phi with p = 1, \beta = 0.2 
%======================================================
% Creating calibration variables
p = 1;
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.2;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = 0.47;
theta_ss_p1 = fsolve(@(theta) jc_theta_ss(theta), theta0);

% Finding the SS wage given theta_ss
w_ss_p1 = w_theta_ss(theta_ss_p1);

% Finding the SS unemployment rate given theta_ss
u_ss_p1 = u_theta_ss(theta_ss_p1);

% Finding the SS job-finding rate given theta_ss
phi_ss_p1 = phi_theta(theta_ss_p1);

%========================================================
% Calculating SS u, w, and \phi with p = 1.1, \beta = 0.2 
%========================================================
% Creating calibration variables
p = 1.1;
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.2;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = 0.47;
theta_ss_p11 = fsolve(@(theta) jc_theta_ss(theta), theta0);

% Finding the SS wage given theta_ss
w_ss_p11 = w_theta_ss(theta_ss_p11);

% Finding the SS unemployment rate given theta_ss
u_ss_p11 = u_theta_ss(theta_ss_p11);

% Finding the SS job-finding rate given theta_ss
phi_ss_p11 = phi_theta(theta_ss_p11);

% To avoid more mindless algebra and derivatives (ASK PROFESSOR MUELLER IF
% THIS IS ACTUALLY THE DESIRED METHOD), we will use the ratio of per cent
% changes definition of elasticity with respect to p.
u_elasticity_wrt_p = ((u_ss_p11 - u_ss_p1)/(u_ss_p11 + u_ss_p1))/((1.1 - 1)/(1.1 + 1));
w_elasticity_wrt_p = ((w_ss_p11 - w_ss_p1)/(w_ss_p11 + w_ss_p1))/((1.1 - 1)/(1.1 + 1));
phi_elasticity_wrt_p = ((phi_ss_p11 - phi_ss_p1)/(phi_ss_p11 + phi_ss_p1))/((1.1 - 1)/(1.1 + 1));
%=======
% ANSWER
%=======
% Using p \in {1, 1.1}, the SS elasticities of unemployment rate, wage,
% and job-finding rate with respect to p are -0.4583, 0.9851, and 0.4708
% respectively.

% We see that the SS unemployment rate in this alternative calibration is
% too low relative to the socially efficient unemployment rate. This is
% because \beta = 0.2 < \eta = 0.72, meaning Hosios' condition is not
% met. Intuitively, this inequality will result in firms creating more
% congestion to other firms compared to congestion created amongst workers
% (not that this is at the margin).
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 2d: Compute the SS elasticities of u, w, and \phi for the same calibration as in 2a., but where you set \beta = 0.2 and where you vary the parameter b in 0.08-increments from 0.40 to 0.96. Plot the SS elasticities in a graph as a function of b.
%=================================================================
% Calculating SS u, w, and \phi with p = 1, \beta = 0.2, varying b 
%=================================================================
% Creating calibration variables
p = 1;
b = (0.40:0.08:0.96);
lambda = 0.03;
eta = 0.72;
beta = 0.2;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = repmat(0.47, 1, length(b));
theta_ss_p1 = fsolve(@(theta) jc_theta_ss(theta), theta0);

% Finding the SS wage given theta_ss
w_ss_p1 = w_theta_ss(theta_ss_p1);

% Finding the SS unemployment rate given theta_ss
u_ss_p1 = u_theta_ss(theta_ss_p1);

% Finding the SS job-finding rate given theta_ss
phi_ss_p1 = phi_theta(theta_ss_p1);

%===================================================================
% Calculating SS u, w, and \phi with p = 1.1, \beta = 0.2, varying b 
%===================================================================
% Creating calibration variables
p = 1.1;
b = (0.40:0.08:0.96);
lambda = 0.03;
eta = 0.72;
beta = 0.2;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = repmat(0.47, 1, length(b));
theta_ss_p11 = fsolve(@(theta) jc_theta_ss(theta), theta0);

% Finding the SS wage given theta_ss
w_ss_p11 = w_theta_ss(theta_ss_p11);

% Finding the SS unemployment rate given theta_ss
u_ss_p11 = u_theta_ss(theta_ss_p11);

% Finding the SS job-finding rate given theta_ss
phi_ss_p1 = phi_theta(theta_ss_p11);

% To avoid more mindless algebra and derivatives (ASK PROFESSOR MUELLER IF
% THIS IS ACTUALLY THE DESIRED METHOD), we will use the ratio of per cent
% changes definition of elasticity with respect to p.
u_elasticity_wrt_p = ((u_ss_p11 - u_ss_p1)./(u_ss_p11 + u_ss_p1))./((1.1 - 1)./(1.1 + 1));
w_elasticity_wrt_p = ((w_ss_p11 - w_ss_p1)./(w_ss_p11 + w_ss_p1))./((1.1 - 1)./(1.1 + 1));
phi_elasticity_wrt_p = ((phi_ss_p11 - phi_ss_p1)./(phi_ss_p11 + phi_ss_p1))./((1.1 - 1)./(1.1 + 1));

% Plotting SS elasticities wrt p
plot(b, u_elasticity_wrt_p, 'LineWidth', 2);
hold on
plot(b, w_elasticity_wrt_p, 'LineWidth', 2);
hold on
plot(b, phi_elasticity_wrt_p, 'LineWidth', 2);
grid on
xlabel('b');
ylabel('Elasticity');
hold on

% Creating legend
legend('SS unemployment rate (Blue)', 'Wage (Red)', 'Job-finding rate (Yellow)', 'Location', 'Northwest');

% Creating title
title({'SS Elasticities with respect to p,', 'unemployment rate, wages, and the job-finding rate', 'b \in \{0.40, 0.48, \ldots 0.96\}'});

saveas(gcf, 'path\to\graphics\2d_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 2e: Assume now that the wage is set according to the following wage norm: w = \alpha*\omega + (1 - \alpha)*w_{NB}, where w_{NB} is the Nash-bargained wage in the standard model in 2a. Calibrate \omega such that \omega = w_{NB} = w in the SS equilibrium for p = 1. Compute and plot the SS elasticities of u, w, and \phi with respect to p as afunction of \alpha (going from 0 to 1 in 0.1 increments).
%======================================================================
% Calculating SS u, w, and \phi with p = 1, \beta = 0.2, varying \alpha 
%======================================================================
% Creating calibration variables
p = 1;
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;
% Setting parameter \alpha
alpha = (0:0.1:1);

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS Nash-bargained wage equation as a function of \theta
w_theta_nb_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

%=====
% NOTE
%=====
% We want to calibrate \omega such that the new wage equation w equals both
% \omega and the Nash-bargained wage w_{NB} in the SS equilibrium with p =
% 1. Recall that in the calibraiton of part 2a, we found that w_{NB} =
% 0.9871 in SS equilibrium. We use this value to calibrate \omega.
%=========
% END NOTE
%=========
% Setting omega to equal w_{NB} in SS equilibrium from part 2a
omega = 0.9871;

% Creating SS wage equation that is a weighted average of \omega and
% w_theta_nb_ss
w_theta_ss = @(theta) alpha.*omega + (1 - alpha).*w_theta_nb_ss(theta);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = repmat(0.47, 1, length(alpha));
theta_ss_p1 = fsolve(@(theta) jc_theta_ss(theta), theta0);

% Finding the SS wage given theta_ss
w_ss_p1 = w_theta_ss(theta_ss_p1);

% Finding the SS unemployment rate given theta_ss
u_ss_p1 = u_theta_ss(theta_ss_p1);

% Finding the SS job-finding rate given theta_ss
phi_ss_p1 = phi_theta(theta_ss_p1);
%=======
% ANSWER
%=======
% When checking, we see that the SS equilibrium wage equals both \omega and
% the SS equilibrium Nash-bargained wage (all three in turn equal to
% 0.9871, which is the equilibrium wage found in part 2a).
%===========
% END ANSWER
%===========

%========================================================================
% Calculating SS u, w, and \phi with p = 1.1, \beta = 0.2, varying \alpha 
%========================================================================
% Creating calibration variables
p = 1.1;
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;
% Setting parameter \alpha
alpha = (0:0.1:1);

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS Nash-bargained wage equation as a function of \theta
w_theta_nb_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

%=====
% NOTE
%=====
% We want to calibrate \omega such that the new wage equation w equals both
% \omega and the Nash-bargained wage w_{NB} in the SS equilibrium with p =
% 1. Recall that in the calibraiton of part 2a, we found that w_{NB} =
% 0.9871 in SS equilibrium. We use this value to calibrate \omega.
%=========
% END NOTE
%=========
% Setting omega to equal w_{NB} in SS equilibrium from part 2a
omega = 0.9871;

% Creating SS wage equation that is a weighted average of \omega and
% w_theta_nb_ss
w_theta_ss = @(theta) alpha.*omega + (1 - alpha).*w_theta_nb_ss(theta);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = repmat(0.47, 1, length(alpha));
theta_ss_p11 = fsolve(@(theta) jc_theta_ss(theta), theta0);

% Finding the SS wage given theta_ss
w_ss_p11 = w_theta_ss(theta_ss_p11);

% Finding the SS unemployment rate given theta_ss
u_ss_p11 = u_theta_ss(theta_ss_p11);

% Finding the SS job-finding rate given theta_ss
phi_ss_p11 = phi_theta(theta_ss_p11);

% To avoid more mindless algebra and derivatives (ASK PROFESSOR MUELLER IF
% THIS IS ACTUALLY THE DESIRED METHOD), we will use the ratio of per cent
% changes definition of elasticity with respect to p.
u_elasticity_wrt_p = ((u_ss_p11 - u_ss_p1)./(u_ss_p11 + u_ss_p1))./((1.1 - 1)./(1.1 + 1));
w_elasticity_wrt_p = ((w_ss_p11 - w_ss_p1)./(w_ss_p11 + w_ss_p1))./((1.1 - 1)./(1.1 + 1));
phi_elasticity_wrt_p = ((phi_ss_p11 - phi_ss_p1)./(phi_ss_p11 + phi_ss_p1))./((1.1 - 1)./(1.1 + 1));

% Plotting SS elasticities wrt p
plot(alpha, u_elasticity_wrt_p, 'LineWidth', 2);
hold on
plot(alpha, w_elasticity_wrt_p, 'LineWidth', 2);
hold on
plot(alpha, phi_elasticity_wrt_p, 'LineWidth', 2);
grid on
xlabel('b');
ylabel('Elasticity');
hold on

% Creating legend
legend('SS unemployment rate (Blue)', 'Wage (Red)', 'Job-finding rate (Yellow)', 'Location', 'Northwest');

% Creating title
title({'SS Elasticities with respect to p,', 'unemployment rate, wages, and the job-finding rate', '\alpha \in \{0, 0.1, \ldots 1\}'});

saveas(gcf, 'path\to\graphics\2e_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 2f: Going back to the model in part 2a, plot the SS relationship between v and u for different values of p in the interval from 0.8 to 1.2. On a different graph, plot the SS relationship betwen v and u for different values of the separation rate \lambda (p = 1). Which plot fits the empirical relationship job openings and unemployment better?
%===========================================================
% Calculating SS u and v with part 2a calibration, varying p
%===========================================================
% Creating calibration variables
p = (0.8:0.05:1.2);
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = repmat(0.47, 1, length(p));
theta_ss = fsolve(@(theta) jc_theta_ss(theta), theta0);
%theta_ss = fsolve(@(theta) ((p - (1 - beta).*b - beta.*p)/beta.*k) - ((r + lambda)./(beta.*q_theta(theta))) - theta, theta0);

% Finding the SS unemployment rate given theta_ss
u_ss = u_theta_ss(theta_ss);

% Finding the SS job vacancy rate given theta_ss
v_ss = theta_ss.*u_ss;

% Plotting Beveridge diagram
plot(u_ss*100, v_ss*100, 'LineWidth', 2);
hold on
grid on
xlabel('Unemployment rate (%)');
ylabel('Job openings rate (%)');
hold on

% Creating title
title({'SS Beveridge curve diagram (u - v plane),', 'p \in \{0.8, 0.85, \ldots 1.2\}'});

saveas(gcf, 'path\to\graphics\2fi_plot.png');
close(gcf);

%=================================================================
% Calculating SS u and v with part 2a calibration, varying \lambda
%=================================================================
% Creating calibration variables
p = 1;
b = 0.4;
lambda = (0.03:0.01:0.27);
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Creating SS wage equation as a function of \theta
w_theta_ss = @(theta) (1 - beta)*b + beta.*(p + theta.*k);

% Creating SS job creation condition equation as a function of \theta
jc_theta_ss = @(theta) ((p - w_theta_ss(theta))./(r + lambda)) - k./q_theta(theta);

% Creating SS Beveridge curve equation as a function of \theta
u_theta_ss = @(theta) lambda./(lambda + phi_theta(theta));

% Finding \theta such that SS job creation condition (with wage equation
% plugged in) is equal to zero
theta0 = repmat(0.47, 1, length(lambda));
theta_ss = fsolve(@(theta) jc_theta_ss(theta), theta0);
%theta_ss = fsolve(@(theta) ((p - (1 - beta).*b - beta.*p)/beta.*k) - ((r + lambda)./(beta.*q_theta(theta))) - theta, theta0);

% Finding the SS unemployment rate given theta_ss
u_ss = u_theta_ss(theta_ss);

% Finding the SS job vacancy rate given theta_ss
v_ss = theta_ss.*u_ss;

% Plotting Beveridge diagram
plot(u_ss*100, v_ss*100, 'LineWidth', 2);
hold on
grid on
xlabel('Unemployment rate (%)');
ylabel('Job openings rate (%)');
hold on

% Creating title
title({'SS Beveridge curve diagram (u - v plane),', '\lambda \in \{0.03, 0.04, \ldots 0.27\}'});

saveas(gcf, 'path\to\graphics\2fii_plot.png');
close(gcf);
%=======
% ANSWER
%=======
% Based on visual comparison, we see that the Beveridge curve diagram from
% varying productivity p results in a plot that fits the empirical
% relationship between job opens and unemployment better.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 2g: Model with endogenous separations. Assume that the matching function takes the form M(u,v) = \mu*(u^{\eta})*(v^{1 - \eta}), the distribution of match-specific productivity draws x follows U(0, 1), and calibrate the model at the monthly level by choosing the following parameter values: p = 1, b = 0.4, \eta = 0.72, \beta = 0.72, and r = 0.004. Set \mu = 1 and set the remaining parameter k such that the monthly job finding rate, \phi(\theta), is 0.30 and lambda such that the monthly separation rate is 0.03 in the SS equilibrium.
% Creating calibration variables
p = 1;
b = 0.4;
% Creating vector of \lambda from 0 to 1 in 0.00001 increments
lambda = (0:0.0001:1);
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*(theta.^(-eta));

% Creating job-finding rate function, \phi(\theta)
phi_theta = @(theta) theta.*q_theta(theta);

% Solving for theta that yields job openings rate of 0.3
theta0 = 0.47;
theta_for_phi03_ss = fsolve(@(theta) phi_theta(theta) - 0.3, theta0);

% Creating k function, which is a function of R. This is created by
% rearranging the job creation condition such that k is on the LHS
k_R_ss = @(R) ((1 - beta)*p*(1 - R)*q_theta(theta_for_phi03_ss))./(r + lambda);

% Creating job destruction condition as a function of R
jd_R_ss = @(R) R - (b/p) - (beta/(1 - beta)).*((theta_for_phi03_ss.*k_R_ss(R))./p) + (lambda./(r + lambda)).*((1 - R)/2);
  
% Setting initial value of R
R0 = repmat(0.47, 1, length(lambda));
R_ss_vector = fsolve(@(R) jd_R_ss(R), R0);

% Creating vector of separation rates \lambda*G(R)
separation_rate_ss = lambda'.*unifcdf(R_ss_vector)';

% Creating matrix that pairs separation rates with corresponding \lambda
% and SS equilibrium R
l_R_sr_matrix = [lambda' R_ss_vector' separation_rate_ss];

% Adding absolute difference between separation rate and 0.03 target as
% leftmost column in the above matrix
l_R_sr_diff_matrix = [l_R_sr_matrix abs(separation_rate_ss - 0.03)];

% Finding the \lambda value and SS equilibrium R that results in a
% separation rate to be closest to the 0.03 target
min_diff = min(l_R_sr_diff_matrix(:, 4));
[row, column] = find(l_R_sr_diff_matrix == min_diff);

% Obtaining value for SS equilibrium R given the separation rate target of
% 0.03.
R_ss = R_ss_vector(row);
%=======
% ANSWER
%=======
% We see that when given \lambda = 0.0331, the reservation productivity R
% in SS equilibrium is 0.9059.
%===========
% END ANSWER
%===========

% Obtaining value for paramter k
k_ss_vector = k_R_ss(R_ss_vector);
k_ss = k_ss_vector(row);
%=======
% ANSWER
%=======
% With our specified calibration, it's found that k = 15.7045.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 2h: For the model calibrated in 2g, compute the SS elasticity of u, w, and \phi with respect to p. Moreover, plot the SS relationship between \theta and u (i.e., the Beveridge curve diagram) for different values of p in the interval from 0.8 to 1.2.
%=====
% NOTE
%=====
% Observe that part 2g was to figure out every calibrated variable for the
% model. Specifically, we solved for parameter k, given the theta value
% that results in \phi = 0.3. To solve for k, we obtained the SS
% equilibrium value of R and parameter value \lambda.

% What this translates for part 2h is that our model with endogenous
% separations is now fully parameterised: We have parameters lambda and k.
% We now need to resolve the model with p = 1 and p = 1.1.
%=========
% END NOTE
%=========
%=========================================
% Calculating SS u, w, and \phi with p = 1
%=========================================
% Creating calibration variables
p = 1;
b = 0.4;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameters found in part 2g
lambda = 0.0331;
k = k_ss;

%=====
% NOTE
%=====
% x(1) = R;
% x(2) = \theta.
%=========
% END NOTE
%=========
% Creating k function, which is a function of R and \theta. This is created
% by rearranging the job creation condition such that k is on the LHS
k_R_theta_ss = @(x) ((1 - beta)*p*(1 - x(1))*q_theta(x(2)))./(r + lambda);

% Creating job destruction condition as a function of R and \theta
jd_R_theta_ss = @(x) x(1) - (b/p) - (beta/(1 - beta)).*((x(2).*k_R_theta_ss(x))./p) + (lambda./(r + lambda)).*((1 - x(1))/2);

% Creating SS equilibrium Beveridge curve as a function of R and \theta
u_R_theta_ss = @(R, theta) lambda*unifcdf(R)/(lambda*unifcdf(R) + phi_theta(theta));

%=====
% NOTE
%=====
% Because we are given that match-specific productivity draws x follows the
% distribution U(0, 1), we will use the expected value of x for determining
% the SS equilibrium elasticity of wage with respect to p. As a result, we
% will define the wage function using E[x] = (1 + R)/2.
%=========
% END NOTE
%=========
w_R_theta_ss = @(R, theta) (1 - beta)*b + beta*(p*((1 + R)/2) + theta*k);

% Setting initial values of R and \lambda
R0 = 0.47;
theta0 = 0.47;
initial = [R0 theta0];
R_theta_vector = fsolve(@(x) jd_R_theta_ss(x), initial);

% Obtaining equilibrium u, R, w, and \phi
u_ss_p1 = u_R_theta_ss(R_theta_vector(1), R_theta_vector(2));
R_ss_p1 = R_theta_vector(1);
w_ss_p1 = w_R_theta_ss(R_theta_vector(1), R_theta_vector(2));
phi_ss_p1 = phi_theta(R_theta_vector(2));

%===========================================
% Calculating SS u, w, and \phi with p = 1.1
%===========================================
% Creating calibration variables
p = 1.1;
b = 0.4;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameters found in part 2g
lambda = 0.0331;
k = k_ss;

%=====
% NOTE
%=====
% x(1) = R;
% x(2) = \theta.
%=========
% END NOTE
%=========
% Creating k function, which is a function of R and \theta. This is created
% by rearranging the job creation condition such that k is on the LHS
k_R_theta_ss = @(x) ((1 - beta)*p*(1 - x(1))*q_theta(x(2)))./(r + lambda);

% Creating job destruction condition as a function of R and \theta
jd_R_theta_ss = @(x) x(1) - (b/p) - (beta/(1 - beta)).*((x(2).*k_R_theta_ss(x))./p) + (lambda./(r + lambda)).*((1 - x(1))/2);

% Creating SS equilibrium Beveridge curve as a function of R and \theta
u_R_theta_ss = @(R, theta) lambda*unifcdf(R)/(lambda*unifcdf(R) + phi_theta(theta));

%=====
% NOTE
%=====
% Because we are given that match-specific productivity draws x follows the
% distribution U(0, 1), we will use the expected value of x for determining
% the SS equilibrium elasticity of wage with respect to p. As a result, we
% will define the wage function using E[x] = (1 + R)/2.
%=========
% END NOTE
%=========
w_R_theta_ss = @(R, theta) (1 - beta)*b + beta*(p*((1 + R)/2) + theta*k);

% Setting initial values of R and \lambda
R0 = 0.47;
theta0 = 0.47;
initial = [R0 theta0];
R_theta_vector = fsolve(@(x) jd_R_theta_ss(x), initial);

% Obtaining equilibrium u, R, w, and \phi
u_ss_p11 = u_R_theta_ss(R_theta_vector(1), R_theta_vector(2));
R_ss_p11 = R_theta_vector(1);
w_ss_p11 = w_R_theta_ss(R_theta_vector(1), R_theta_vector(2));
phi_ss_p11 = phi_theta(R_theta_vector(2));

% To avoid more mindless algebra and derivatives (ASK PROFESSOR MUELLER IF
% THIS IS ACTUALLY THE DESIRED METHOD), we will use the ratio of per cent
% changes definition of elasticity with respect to p.
u_elasticity_wrt_p = ((u_ss_p11 - u_ss_p1)./(u_ss_p11 + u_ss_p1))./((1.1 - 1)./(1.1 + 1));
w_elasticity_wrt_p = ((w_ss_p11 - w_ss_p1)./(w_ss_p11 + w_ss_p1))./((1.1 - 1)./(1.1 + 1));
phi_elasticity_wrt_p = ((phi_ss_p11 - phi_ss_p1)./(phi_ss_p11 + phi_ss_p1))./((1.1 - 1)./(1.1 + 1));
%=======
% ANSWER
%=======
% Using p \in {1, 1.1}, the SS elasticities of unemployment rate, wage,
% and job-finding rate with respect to p are -0.0307, 0.1746, and 0.0053
% respectively.
%===========
% END ANSWER
%===========

%========================================
% Calculating SS u, and \theta, varying p
%========================================
% Creating calibration variables
p = (0.8:0.05:1.2);
b = 0.4;
eta = 0.72;
beta = 0.72;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameters found in part 2g
lambda = 0.0331;
k = k_ss;

% Looping and solving the model for each value of p
for i = 1:length(p)
  %=====
  % NOTE
  %=====
  % x(1) = R;
  % x(2) = \theta.
  %=========
  % END NOTE
  %=========
  % Creating k function, which is a function of R and \theta. This is created
  % by rearranging the job creation condition such that k is on the LHS
  k_R_theta_ss = @(x) ((1 - beta)*p(i)*(1 - x(1))*q_theta(x(2)))./(r + lambda);

  % Creating job destruction condition as a function of R and \theta
  jd_R_theta_ss = @(x) x(1) - (b/p(i)) - (beta/(1 - beta)).*((x(2).*k_R_theta_ss(x))./p(i)) + (lambda./(r + lambda)).*((1 - x(1))/2);

  % Creating SS equilibrium Beveridge curve as a function of R and \theta
  u_R_theta_ss = @(R, theta) lambda*unifcdf(R)/(lambda*unifcdf(R) + phi_theta(theta));

  % Setting initial values of R and \lambda
  R0 = 0.47;
  theta0 = 0.47;
  initial = [R0 theta0];
  R_theta_vector = fsolve(@(x) jd_R_theta_ss(x), initial);

  % Storing SS equilibrium unemployment rate and \theta
  u_theta_ss_matrix(i, 1) = u_R_theta_ss(R_theta_vector(1), R_theta_vector(2));
  u_theta_ss_matrix(i, 2) = R_theta_vector(2);
end

% Plotting Beveridge diagram
plot(u_theta_ss_matrix(1:length(p), 2), u_theta_ss_matrix(1:length(p), 1)*100, 'LineWidth', 2);
hold on
grid on
xlabel('Labour market tightness');
ylabel('Unemployment rate (%)');
hold on

% Creating title
title({'SS Beveridge curve diagram (\theta - u plane),', 'p \in \{0.8, 0.85, \ldots 1.2\}'});

saveas(gcf, 'path\to\graphics\2h_plot.png');
close(gcf);