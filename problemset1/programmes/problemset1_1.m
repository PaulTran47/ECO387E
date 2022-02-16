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
% ECO387E Problem Set 1, 1
% Paul Le Tran, plt377
% 13 February, 2022
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
%% Part 1e: Plot the density f(w) and the values \overline{w} and w_{R} for the following calibration of the model: p = 1, b - 0.7, \lambda_{u} = 0.4, \lambda_{e} = 0.1, \delta = 0.03
% Creating calibration variables
p_1e = 1;
b = 0.7;
lambda_u = 0.4;
lambda_e = 0.1;
delta = 0.03;
k_u = lambda_u/delta;
k_e = lambda_e/delta;

% Creating w_{R} with calibration
w_r_1e = (((1 + k_e)^2)*b + (k_u - k_e)*k_e*p_1e)/(((1 + k_e)^2) + (k_u - k_e)*k_e);

% Creating \overline{w} with calibration
w_bar_1e = p_1e - (p_1e - w_r_1e)/((1 + k_e)^2);

% Creating density function f(w)
f_w_1e = @(w) ((1 + k_e)/(2*k_e))/sqrt((p_1e - w_r_1e)*(p_1e - w));

% Plotting f(w), w_bar, and w_r
fplot(f_w_1e, [w_r_1e - 0.05, w_bar_1e + 0.05], 'LineWidth', 2);
hold on
xline(w_r_1e, 'LineWidth', 1.5);
xline(w_bar_1e, 'LineWidth', 1.5);
grid on
xlabel('Wage w');
ylabel('Density');
hold on
title({'Density function f(w),', 'w \in [0.8919, 0.9942]'});

saveas(gcf, 'path\to\graphics\1e_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 1f: Plot density in 1e as well as the density f(w) and the values of \overline{w} and w_{R} for the following calibration: p = 1.1, b = 0.7, \lambda_{u} = 0.4, \lambda_{e} = 0.1, \delta = 0.03. Describe in words how a higher level of aggregate productivity changes the wage distribution.
% Observe that the only change in calibration variables is p
p_1f = 1.1;

% Creating w_{R} with calibration
w_r_1f = (((1 + k_e)^2)*b + (k_u - k_e)*k_e*p_1f)/(((1 + k_e)^2) + (k_u - k_e)*k_e);

% Creating \overline{w} with calibration
w_bar_1f = p_1f - (p_1f - w_r_1f)/((1 + k_e)^2);

% Creating density function f(w)
f_w_1f = @(w) ((1 + k_e)/(2*k_e))/sqrt((p_1f - w_r_1e)*(p_1f - w));

% Plotting f(w) for 1e, f(w) for 1f, and w_bar and w_r for 1f
fplot(f_w_1e, [w_r_1e - 0.05, w_bar_1e + 0.05], 'LineWidth', 2);
hold on
fplot(f_w_1f, [w_r_1f - 0.05, w_bar_1f + 0.05], 'LineWidth', 2);
hold on
xline(w_r_1f, 'LineWidth', 1.5);
xline(w_bar_1f, 'LineWidth', 1.5);
grid on
xlabel('Wage w');
ylabel('Density');
hold on
title({'Density function f(w),', 'w \in [0.8919, 0.9942] for 1e (not shown)', 'w \in [0.9559, 1.0923] for 1f'});

saveas(gcf, 'path\to\graphics\1f_plot.png');
close(gcf);
%=======
% ANSWER
%=======
% From our plot comparing the two density functions, we see that a higher
% level of aggregate productivity changes the wage distribution by shifting
% it foward. This change also occurs to both the reservation wage and
% maximum possible wage.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 1g: Plot the reservation wage, w_{R}, for different values of \lambda_{e} in the range from 0 to 0.4 (otherwqise follow the same calibration as in 1e). What is the main intuition for why the reservation wage is decreasing in \lambda_{e}?
% Creating calibration variables
p = 1;
b = 0.7;
lambda_u = 0.4;
lambda_e = 0.1;
delta = 0.03;
k_u = lambda_u/delta;
% Because \lambda_e is now a changing variable, we will make k_e a function
k_e = @(lambda) lambda/delta;

% Creating w_{R} function to take values of \lambda_{e}
w_r = @(lambda) (((1 + k_e(lambda))^2)*b + (k_u - k_e(lambda))*k_e(lambda)*p)/(((1 + k_e(lambda))^2) + (k_u - k_e(lambda))*k_e(lambda));

% Plotting reservation wage for \lambda_{e} \in [0, 0.4]
fplot(w_r, [0, 0.4], 'LineWidth', 2);
grid on
xlabel('\lambda_{e}');
ylabel('Wage');
hold on
title({'Reservation wage w_{R} as a function of \lambda_{e},', '\lambda_{e} \in [0, 0.4]'});

saveas(gcf, 'path\to\graphics\1g_plot.png');
close(gcf);
%=======
% ANSWER
%=======
% From our plot, we see that a employed-person job offer arrival rate
% of zero will result in the employed person to prefer receiving
% unemployment benefits b as a minimum. This is because with the small
% chance of their job being destroyed, they would just settle on setting
% their reservation wage to be $b$ (as an unemployed person would do in the
% model).

% However, as the job arrival rate for employed people being to increase,
% early positive values will result in an employed person to quickly
% increase their reservation wage. This is because as the rate increases,
% they know their odds of receiving a job offer with better wage increases.
% Therefore, they raise their minimum threshold/standard for accepting a
% job offer.

% The reservation wage will reach a maximum as $\lambda_{e}$ continues to
% increase, however. Once so, the reservation wage begins to decrease back
% to the value of b at $\lambda_{e} = 0.4$. Note that this decrease will
% continue as the job arrival rate for employed people increase.
% Intuitively, this trend is because NEED TO FIGURE OUT WHY.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 1i: Plot the density f(w(p)) for the following calibration of the model: p_bar = 1, b = 0.7, \lambda_{u} = \lambda_{e} = 0.4, \delta = 0.03. On the same graph, plot the density f(w) of the model without heterogeneous firms, for the following calibration: p = 1, b= 0.7, \lambda_{u} = \lambda_{e} = 0.4, \delta = 0.03. In words, briefly describe the differences in the distributions shown.
% Creating calibration variables for f(w(p))
p_bar = 1;
b = 0.7;
lambda = 0.4;
delta = 0.03;
k = lambda/delta;

% Creating wage function w(p) and derivative w'(p)
w_p = @(p) p - ((p_bar - b + k*p_bar - k*p)/k)*(1 - ((p_bar - b + k*p_bar - k*p)/(p_bar - b + k*p_bar - b*p)));
w1_p = @(p) (p - w_p(p))*((2*k*unifpdf(p, b, p_bar))/(1 + k*(1 - unifcdf(p, b, p_bar))));

% Creating density function in the model with heterogeneous firms f(w(p))
f_w_p_1i = @(p) unifpdf(p, b, p_bar)/w1_p(p);

% Creating calibration variables for f(w)
p = 1;
b = 0.7;
lambda = 0.4;
delta = 0.03;
k = lambda/delta;

% Creating w_{R} with calibration
w_r = (((1 + k)^2)*b)/(((1 + k)^2));

% Creating \overline{w} with calibration
w_bar = p - (p - w_r)/((1 + k)^2);

% Creating density function f(w)
f_w = @(w) ((1 + k)/(2*k))/sqrt((p - w_r)*(p - w));

% Plotting f(w(p))
fplot(f_w_p_1i, [b, p_bar], 'LineWidth', 2);
grid on
xlabel('Productivity p');
ylabel('Density');
hold on
title({'Density function f(w(p)) (model with heterogenous firms),', 'p \in [0.7, 1]'});

saveas(gcf, 'path\to\graphics\1ia_plot.png');
close(gcf);

% Plotting f(w)
fplot(f_w, [b, p - 0.0014], 'LineWidth', 2);
grid on
xlabel('Wage w');
ylabel('Density');
hold on
title({'Density function f(w) (model with homogenous firms),', 'w \in [0.7, 1]'});

saveas(gcf, 'path\to\graphics\1ib_plot.png');
close(gcf);
%=======
% ANSWER
%=======
% When comparing the two density functions together, we initially see that
% the curve/shape of both distributions differ greatly. For $f(w(p))$ in
% the model with heterogeneous firms, the density function is slightly
% convex upward, whilst the function is also slightly decreasing. For f(w)
% in the model with homogeneous firms, the density function takes the shape
% as ones found in previous parts: Convex as well but with an increasing
% exponential shape.

% Furthermore, the range of densities that each function can take differ as
% well. $f(w(p))$ for the given calibration is in the single digits, whilst
% $f(w)$ exponentially grows towards infinity as wage increases.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 1h: Model with heterogeneous firms with J(p) following the truncated normal distribution with lower bound b = 0.7, upper bound p_bar = 1, mean \mu = 0.85, and standard deviation \sigma = 0.15. Compute w(p) for 31 points on the interval between 0.7 and 1 for the same calibration as in 1i. Plot the density f(w(p)) for each of these 31 points. On the same graph, plot the density for the model with heterogeneous against solved in 1i. In words, briefly describe the differences in the distributions.
% Creating calibration variables for f(w(p))
p_bar = 1;
b = 0.7;
lambda = 0.4;
delta = 0.03;
k = lambda/delta;

% Creating a function for the truncated normal distribution with given
% specifications.
% Setting pre-truncated normal distribution parameters
pretrunc_mu = 0.85;
pretrunc_sigma = 0.15;

% Creating untruncated normal distribution
untruncated = makedist('Normal', pretrunc_mu, pretrunc_sigma);

% Creating J(p)
J_p = truncate(untruncated, pretrunc_mu - 0.15, pretrunc_mu + 0.15);

% Creating wage function w(p) and derivative w'(p)
% Creating integrand
integrand = @(p) (1./(1 + k.*(1 - cdf(J_p, p)))).^2;
% Creating integral inside of w_p
w_p_integral = @(p) integral(integrand, 0.7, p);
w_p = @(p) p - ((1 + k.*(1 - cdf(J_p, p))).^2).*w_p_integral(p);
w1_p = @(p) (p - w_p(p))*((2*k*pdf(J_p, p))/(1 + k*(1 - cdf(J_p, p))));

% Creating density function in the model with heterogeneous firms f(w(p))
f_w_p_1h = @(p) pdf(J_p, p)./w1_p(p);

% Plotting f(w(p)) for uniform and truncated normal distributions
fplot(f_w_p_1h, [b + 0.01, p_bar], 'LineWidth', 2);
hold on
fplot(f_w_p_1i, [b + 0.01, p_bar], 'LineWidth', 2);
grid on
xlabel('Productivity p');
ylabel('Density');
hold on

% Creating legend
legend('Trunctated normal (Blue)', 'Uniform (Red)', 'Location', 'Northeast');

% Creating title
title({'Density function f(w(p)) (model with heterogenous firms),', 'p \in \{0.71, 0.72, \ldots 1\}'});

saveas(gcf, 'path\to\graphics\1h_plot.png');
close(gcf);
%=======
% ANSWER
%=======
% When comparing densities $f(w(p))$ from the model with heterogeneous
% firms but with $J(p)$ following a uniform and truncated normal
% distributions (both having the same calibration), we see immediately that
% the density using the truncated normal distribution is much more convex.
% We also see that as productivity approaches values close to
% $\overline{p} = 1$, the two density functions converge together. However,
% productivity values close to $b = 0.7$ see the density using the
% truncated normal distribution to be much greater than that using the
% uniform distribution.
%===========
% END ANSWER
%===========