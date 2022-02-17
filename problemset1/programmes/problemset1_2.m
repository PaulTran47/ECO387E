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
q_theta = @(theta) mu.*((1./theta).^eta);

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
% Given our calibration values and mismatching parameter $\mu = 1$, we find
% that $k = 1.6043% will result in SS equilibrium $u = 0.05$, $\theta = 0.1343$ w = 0.9871$, $\phi(\theta) = 0.57$.
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
q_theta = @(theta) mu.*((1./theta).^eta);

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
% With productivity $p = 1.1$ and every other calibration parameter set the
% same as in part 2a., we find that SS $u = 0.0479, \phi(\theta) = 0.5957, w = 1.0856$.
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
% Using $p \in {1, 1.1}$, the SS elasticities of unemployment rate, wage,
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
q_theta = @(theta) mu.*((1./theta).^eta);

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
q_theta = @(theta) mu.*((1./theta).^eta);

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
% Using $p \in {1, 1.1}$, the SS elasticities of unemployment rate, wage,
% and job-finding rate with respect to p are -0.4583, 0.9851, and 0.4708
% respectively.

% We see that the SS unemployment rate in this alternative calibration is
% too low relative to the socially efficient unemployment rate. This is
% because $\beta = 0.2 < \eta = 0.72$, meaning Hosios' condition is not
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
b = [0.40:0.08:0.96];
lambda = 0.03;
eta = 0.72;
beta = 0.2;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*((1./theta).^eta);

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
b = [0.40:0.08:0.96];
lambda = 0.03;
eta = 0.72;
beta = 0.2;
r = 0.004;
% Setting matching mismatch parameter, \mu to problem specification
mu = 1;
% Setting parameter k to equal the solved value from part 2a.
k = k_ss_a;

% Creating job-filling rate function, q(\theta)
q_theta = @(theta) mu.*((1./theta).^eta);

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