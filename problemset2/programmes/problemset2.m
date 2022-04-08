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
% ECO387E Problem Set 2
% Paul Le Tran, plt377
% 7 April, 2022
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
%% Part a (setup): Setting up parameters and transition probabilities with provided information
% Setting up calibration variables to be the same as in problem set 1, 2.
% Assuming wages are completely sticky and follow w(p_{i}) = \omega, which
% is the Nash-Bargained (NB) wage in an economy without productivity shocks
omega = 0.9871;
k = 1.6043;
b = 0.4;
lambda = 0.03;
eta = 0.72;
beta = 0.72;
r = 0.004;
pl = 0.98;
ph = 1.02;
% Given that P(l|h) = P(h|l) = 1 - \pi = 0.04). Therefore, we know P(h|h) =
% P(l|l) = \pi = 0.96
pi = 0.96;
% Setting up error tolerance for value function iteration
eps = 1e-8;
%==========================================================================

%==========================================================================
%% Part ai: Picking an initial guess for U_{i}, W_{i}, and J_{i} for i = {l, h}.
Ul = 0;
Uh = 0;
Wl = 0;
Wh = 0;
Vl = 0;
Vh = 0;
Jl = 0;
Jh = 0;
%==========================================================================

%==========================================================================
%% Part aii: Given your guess for J_{i}, impose the zero-profit condition from vacancy creation and solve for \theta(p_{i}) for both states
%=======
% ANSWER
%=======
% Because of free entry, we know by the zero-profit condition that in the
% basic model, V = 0. In our model, this will also mean V_{l} = V_{h} = 0.
% Therefore, we can say that equation 3 becomes

% 0 = -k + \frac{q(\theta_{i})}{1 + r}[\pi J_{i} + (1 - \pi)J_{j}], for i, j \in {h, l};

% Some algebra will result in the following equation:

% q(\theta_{i}) = \frac{(1 + r)k}{\pi J_{i} + (1 - \pi)J_{j}}, for i, j \in {h, l}.

% Recall that q(\theta) = m(\frac{1}{\theta}, 1). Assuming that the
% matching function follows a Cobb-Douglas structure with \mu = 1, we have:

% q(\theta_{i}) = (\frac{1}{\theta_{i}})^{1/\eta} \implies \theta_{i} =
% (\frac{1}{q(\theta_{i})})^{1/\eta}.

% \therefore \theta_{i} = [\frac{\pi J_{i} + (1 - \pi)J_{j}}{(1 + r)k}]^{1/\eta}.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Parts aiii, aiv, and av: Performing value function iteration
% Setting up variables to count the number of iterations performed and
% initial error term
iter = 0;
delta = 10;

while delta > eps
  iter = iter + 1;
    
  Jh_new = ph - omega + ((1 - lambda)/(1 + r))*(pi*Jh + (1 - pi)*Jl); 
  Jl_new = pl - omega + ((1 - lambda)/(1 + r))*(pi*Jl + (1 - pi)*Jh); 

  % If the value function iteration process results in market tightness
  % (in either state) to be negative, we force market tightness to equal
  % zero instead.
  xh = ((pi*Jh_new + (1 - pi)*Jl_new)/((1 + r)*k))^(1/eta);
  if xh < 0
    thetah = 0;
  else 
    thetah = xh;
  end
    
  xl = ((pi*Jl_new + (1 - pi)*Jh_new)/((1 + r)*k))^(1/eta);
  if xl < 0
    thetal = 0;
  else 
    thetal = xl;
  end

  Uh_new = b + ((1 - thetah^(1 - eta))/(1+r))*(pi*Uh+(1 - pi)*Ul) + (thetah^(1-eta)/(1+r))*(pi*Wh + (1 - pi)*Wl);
  Ul_new = b + ((1 - thetal^(1 - eta))/(1+r))*(pi*Ul+(1 - pi)*Uh) + (thetal^(1-eta)/(1+r))*(pi*Wl + (1 - pi)*Wh);
    
  Wh_new = omega + ((1 - lambda)/(1+r))*(pi*Wh + (1 - pi)*Wl) + (lambda/(1 + r))*(pi*Uh + (1 - pi)*Ul);
  Wl_new = omega + ((1 - lambda)/(1+r))*(pi*Wl + (1 - pi)*Wh) + (lambda/(1 + r))*(pi*Ul + (1 - pi)*Uh);
    
  % Updating error term
  delta = sqrt((Jh_new - Jh)^2 + (Jl_new - Jl)^2) + sqrt((Wh_new - Wh)^2+(Wl_new - Wl)^2) + sqrt((Uh_new - Uh)^2+(Ul_new - Ul)^2);
    
  Jh = Jh_new;
  Jl = Jl_new;
  Wh = Wh_new;
  Wl = Wl_new;
  Uh = Uh_new;
  Ul = Ul_new;
    
  disp(iter);
  disp(delta);
end

phih = thetah^(1 - eta);
phil = thetal^(1 - eta);
%=======
% ANSWER
%=======
% After performing value function iteration, we obtain the following:
% U_{l} = 238.8522;
% U_{h} = 239.2865;
% W_{l} = 240.0383;
% W_{h} = 240.1457;
% J_{l} = 0.2010;
% J_{h} = 0.5609;
% \theta_{l} = 0.0612;
% \theta_{h} = 0.2228;
% \phi_{l}(\theta_{l}) = 0.4573;
% \phi_{h}(\theta_{h}) = 0.6568.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part bi: Now simulate the economy for T = 10000 periods, starting with aggregate state i = l, by following the following procedure: Define vectors S, P, \Theta, \Phi, V, and U of length T.
T = 10000;
S = ones(T,1);
P = zeros(T,1);
Theta = zeros(T,1);
Phi = zeros(T,1);
V = zeros(T,1);
U = zeros(T,1);
%==========================================================================

%==========================================================================
%% Part bii: Draw a vector of T random uniform numbers, R. Simulate the aggregate state i = l for T periods and record it in vector S with value 1 if low state or value 2 if high state.
rng('default');
rng(1);
R = rand(T,1);

for i = 2:T
  if R(i) < 1 - pi && S(i - 1) == 1
    S(i) = 2;
  else 
    if R(i) < 1 - pi && S(i - 1) == 2
      S(i) = 1;
    end 
  end
    
  if R(i) >= 1 - pi
    S(i) = S(i - 1);
  end
end
%==========================================================================

%==========================================================================
%% Part biii: Given the time series for the aggregate state, S, fill in the values of p_{i}, \theta_{i}, and \phi(\theta_{i}) in the corresponding vectors P, \Theta, and \Phi for each period t = 1, \ldots, T.
P(S == 1) = pl;
P(S == 2) = ph;

Theta(S == 1) = thetal;
Theta(S == 2) = thetah;

Phi(S == 1) = phil;
Phi(S == 2) = phih;
%==========================================================================

%==========================================================================
%% Part biv: Simulate the time series of the unemployment rate, U, for t > 1 up to period T, by assuming that U(1) = 0.05 and U(t + 1) = U(t) + \lambda(1 - U(t)) - \Phi(t)*U(t). Compute the time series of vacancies, V.
U(1) = 0.05;
V(1) = U(1)*Theta(1);

for i = 2:T
  U(i) = U(i - 1)+lambda*(1 - U(i - 1)) - Phi(i - 1)*U(i - 1);
  V(i) = U(i)*Theta(i);
end
%==========================================================================

%==========================================================================
%% Part bv: Drop the first 1000 observations in each time series and compute the average of P, \Theta, \Phi, V, and U and the standard deviation of the log of P, \Theta, \Phi, V, and U.
P_m = mean(P(1001:end));
P_sd = std(log(P(1001:end)));

Theta_m = mean(Theta(1001:end));
Theta_sd = std(log(Theta(1001:end)));

Phi_m = mean(Phi(1001:end));
Phi_sd = std(log(Phi(1001:end)));

V_m = mean(V(1001:end));
V_sd = std(log(V(1001:end)));

U_m = mean(U(1001:end));
U_sd = std(log(U(1001:end)));
%==========================================================================

%==========================================================================
%% Part c: Show the sample path for \Phi and U for the last 200 observations of the time series. Does the shape of the two time series differ? If so, explain why.
x = linspace(1, 200, 200);
plot(x', Phi((T - 199):T), 'LineWidth', 2);
hold on
plot(x', U((T - 199):T), 'LineWidth', 2);
grid on
xlabel('Last 200 periods');
hold on

% Creating legend
legend('Job-finding rate (Blue)', 'Unemployment rate (Red)', 'Location', 'Northeast');

% Creating title
title({'Time series of job-finding and unemployment rates,', 'for the last 200 observations of the time period.'});

saveas(gcf, 'path\to\graphics\c_plot.png');
close(gcf);
%=======
% ANSWER
%=======
% The two time series indeed differ in shape. Specifically, we see that in
% periods where the unemployment rate is elevated, the job-finding rate in
% the same time periods decreases to a lower level and remains there. This
% makes sense because as the unemployment rate increases, the job-finding
% rate should decrease as less people are able to move out of unemployment.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part d: Set the initial guess in part a to 100 for all value functions and redo the iterative procedure in part a. Do the value functions converge?
Ul_d = 100;
Uh_d = 100;
Wl_d = 100;
Wh_d = 100;
Vl_d = 100;
Vh_d = 100;
Jl_d = 100;
Jh_d = 100;

% Setting up variables to count the number of iterations performed and
% initial error term
iter = 0;
delta = 10;

while delta > eps
  iter = iter + 1;
    
  Jh_d_new = ph - omega + ((1 - lambda)/(1 + r))*(pi*Jh_d + (1 - pi)*Jl_d); 
  Jl_d_new = pl - omega + ((1 - lambda)/(1 + r))*(pi*Jl_d + (1 - pi)*Jh_d); 

  % If the value function iteration process results in market tightness
  % (in either state) to be negative, we force market tightness to equal
  % zero instead.
  xh_d = ((pi*Jh_d_new + (1 - pi)*Jl_d_new)/((1 + r)*k))^(1/eta);
  if xh_d < 0
    thetah_d = 0;
  else 
    thetah_d = xh_d;
  end
    
  xl_d = ((pi*Jl_d_new + (1 - pi)*Jh_d_new)/((1 + r)*k))^(1/eta);
  if xl_d < 0
    thetal_d = 0;
  else 
    thetal_d = xl_d;
  end

  Uh_d_new = b + ((1 - thetah_d^(1 - eta))/(1+r))*(pi*Uh_d+(1 - pi)*Ul_d) + (thetah_d^(1-eta)/(1+r))*(pi*Wh_d + (1 - pi)*Wl_d);
  Ul_d_new = b + ((1 - thetal_d^(1 - eta))/(1+r))*(pi*Ul_d+(1 - pi)*Uh_d) + (thetal_d^(1-eta)/(1+r))*(pi*Wl_d + (1 - pi)*Wh_d);
    
  Wh_d_new = omega + ((1 - lambda)/(1+r))*(pi*Wh_d + (1 - pi)*Wl_d) + (lambda/(1 + r))*(pi*Uh_d + (1 - pi)*Ul_d);
  Wl_d_new = omega + ((1 - lambda)/(1+r))*(pi*Wl_d + (1 - pi)*Wh_d) + (lambda/(1 + r))*(pi*Ul_d + (1 - pi)*Uh_d);
    
  % Updating error term
  delta = sqrt((Jh_d_new - Jh_d)^2 + (Jl_d_new - Jl_d)^2) + sqrt((Wh_d_new - Wh_d)^2+(Wl_d_new - Wl_d)^2) + sqrt((Uh_d_new - Uh_d)^2+(Ul_d_new - Ul_d)^2);
    
  Jh_d = Jh_d_new;
  Jl_d = Jl_d_new;
  Wh_d = Wh_d_new;
  Wl_d = Wl_d_new;
  Uh_d = Uh_d_new;
  Ul_d = Ul_d_new;
    
  disp(iter);
  disp(delta);
end

phih_d = thetah_d^(1 - eta);
phil_d = thetal_d^(1 - eta);
%=======
% ANSWER
%=======
% When comparing the values from this iterative process with the one
% performed in part a, we see that there is no difference between the two
% part results. Therefore, we conclude that the value functions converge.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part e: Set r = 0 and redo the iterative procedure in part a. Do the value functions converge?
r = 0;

Ul_e = 0;
Uh_e = 0;
Wl_e = 0;
Wh_e = 0;
Vl_e = 0;
Vh_e = 0;
Jl_e = 0;
Jh_e = 0;

% Setting up variables to count the number of iterations performed and
% initial error term
iter = 0;
delta = 10;

while delta > eps
  iter = iter + 1;
    
  Jh_e_new = ph - omega + ((1 - lambda)/(1 + r))*(pi*Jh_e + (1 - pi)*Jl_e); 
  Jl_e_new = pl - omega + ((1 - lambda)/(1 + r))*(pi*Jl_e + (1 - pi)*Jh_e); 

  % If the value function iteration process results in market tightness
  % (in either state) to be negative, we force market tightness to equal
  % zero instead.
  xh_e = ((pi*Jh_e_new + (1 - pi)*Jl_e_new)/((1 + r)*k))^(1/eta);
  if xh_e < 0
    thetah_e = 0;
  else 
    thetah_e = xh_e;
  end
    
  xl_e = ((pi*Jl_e_new + (1 - pi)*Jh_e_new)/((1 + r)*k))^(1/eta);
  if xl_e < 0
    thetal_e = 0;
  else 
    thetal_e = xl_e;
  end

  Uh_e_new = b + ((1 - thetah_e^(1 - eta))/(1+r))*(pi*Uh_e+(1 - pi)*Ul_e) + (thetah_e^(1-eta)/(1+r))*(pi*Wh_e + (1 - pi)*Wl_e);
  Ul_e_new = b + ((1 - thetal_e^(1 - eta))/(1+r))*(pi*Ul_e+(1 - pi)*Uh_e) + (thetal_e^(1-eta)/(1+r))*(pi*Wl_e + (1 - pi)*Wh_e);
    
  Wh_e_new = omega + ((1 - lambda)/(1+r))*(pi*Wh_e + (1 - pi)*Wl_e) + (lambda/(1 + r))*(pi*Uh_e + (1 - pi)*Ul_e);
  Wl_e_new = omega + ((1 - lambda)/(1+r))*(pi*Wl_e + (1 - pi)*Wh_e) + (lambda/(1 + r))*(pi*Ul_e + (1 - pi)*Uh_e);
    
  % Updating error term
  delta = sqrt((Jh_e_new - Jh_e)^2 + (Jl_e_new - Jl_e)^2) + sqrt((Wh_e_new - Wh_e)^2+(Wl_e_new - Wl_e)^2) + sqrt((Uh_e_new - Uh_e)^2+(Ul_e_new - Ul_e)^2);
    
  Jh_e = Jh_e_new;
  Jl_e = Jl_e_new;
  Wh_e = Wh_e_new;
  Wl_e = Wl_e_new;
  Uh_e = Uh_e_new;
  Ul_e = Ul_e_new;
    
  disp(iter);
  disp(delta);
end

phih_e = thetah_e^(1 - eta);
phil_e = thetal_e^(1 - eta);
%=======
% ANSWER
%=======
% No, the value functions do not converge.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part f: Set r = 0.004 as before, but set 1 - \pi = 0.2 and redo the iterative procedures in parts a and b. Show the sample path for \Phi and U for the last 200 observations of the time series. Then compute the average of P, \Theta, \Phi, V, and U and the standard deviation of the log of P, \Theta, \Phi, V, and U. Show averages and standard deviations in the sample table as those in part a. How does switching probability 1 - \pi affect the volatility of job finding and unemployment? Give some intuition for the results.
%===============
% Redoing part a
%===============
r = 0.004;
pi = 0.8;

Ul_f = 0;
Uh_f = 0;
Wl_f = 0;
Wh_f = 0;
Vl_f = 0;
Vh_f = 0;
Jl_f = 0;
Jh_f = 0;

% Setting up variables to count the number of iterations performed and
% initial error term
iter = 0;
delta = 10;

while delta > eps
  iter = iter + 1;
    
  Jh_f_new = ph - omega + ((1 - lambda)/(1 + r))*(pi*Jh_f + (1 - pi)*Jl_f); 
  Jl_f_new = pl - omega + ((1 - lambda)/(1 + r))*(pi*Jl_f + (1 - pi)*Jh_f); 

  % If the value function iteration process results in market tightness
  % (in either state) to be negative, we force market tightness to equal
  % zero instead.
  xh_f = ((pi*Jh_f_new + (1 - pi)*Jl_f_new)/((1 + r)*k))^(1/eta);
  if xh_f < 0
    thetah_f = 0;
  else 
    thetah_f = xh_f;
  end
    
  xl_f = ((pi*Jl_f_new + (1 - pi)*Jh_f_new)/((1 + r)*k))^(1/eta);
  if xl_f < 0
    thetal_f = 0;
  else 
    thetal_f = xl_f;
  end

  Uh_f_new = b + ((1 - thetah_f^(1 - eta))/(1+r))*(pi*Uh_f+(1 - pi)*Ul_f) + (thetah_f^(1-eta)/(1+r))*(pi*Wh_f + (1 - pi)*Wl_f);
  Ul_f_new = b + ((1 - thetal_f^(1 - eta))/(1+r))*(pi*Ul_f+(1 - pi)*Uh_f) + (thetal_f^(1-eta)/(1+r))*(pi*Wl_f + (1 - pi)*Wh_f);
    
  Wh_f_new = omega + ((1 - lambda)/(1+r))*(pi*Wh_f + (1 - pi)*Wl_f) + (lambda/(1 + r))*(pi*Uh_f + (1 - pi)*Ul_f);
  Wl_f_new = omega + ((1 - lambda)/(1+r))*(pi*Wl_f + (1 - pi)*Wh_f) + (lambda/(1 + r))*(pi*Ul_f + (1 - pi)*Uh_f);
    
  % Updating error term
  delta = sqrt((Jh_f_new - Jh_f)^2 + (Jl_f_new - Jl_f)^2) + sqrt((Wh_f_new - Wh_f)^2+(Wl_f_new - Wl_f)^2) + sqrt((Uh_f_new - Uh_f)^2+(Ul_f_new - Ul_f)^2);
    
  Jh_f = Jh_f_new;
  Jl_f = Jl_f_new;
  Wh_f = Wh_f_new;
  Wl_f = Wl_f_new;
  Uh_f = Uh_f_new;
  Ul_f = Ul_f_new;
    
  disp(iter);
  disp(delta);
end

phih_f = thetah_f^(1 - eta);
phil_f = thetal_f^(1 - eta);
%=======
% ANSWER
%=======
% After performing value function iteration, we obtain the following:
% U_{l} = 239.1587;
% U_{h} = 239.2232;
% W_{l} = 240.3013;
% W_{h} = 240.3040;
% J_{l} = 0.2300;
% J_{h} = 0.3252;
% \theta_{l} = 0.0725;
% \theta_{h} = 0.0966;
% \phi_{l}(\theta_{l}) = 0.4797;
% \phi_{h}(\theta_{h}) = 0.5198.
%===========
% END ANSWER
%===========

%======================
% Redoing parts b and c
%======================
S_f = ones(T,1);
P_f = zeros(T,1);
Theta_f = zeros(T,1);
Phi_f = zeros(T,1);
V_f = zeros(T,1);
U_f = zeros(T,1);

rng('default');
rng(1);
R_f = rand(T,1);

for i = 2:T
  if R_f(i) < 1 - pi && S_f(i - 1) == 1
    S_f(i) = 2;
  else 
    if R_f(i) < 1 - pi && S_f(i - 1) == 2
      S_f(i) = 1;
    end 
  end
    
  if R_f(i) >= 1 - pi
    S_f(i) = S_f(i - 1);
  end
end

P_f(S_f == 1) = pl;
P_f(S_f == 2) = ph;

Theta_f(S_f == 1) = thetal_f;
Theta_f(S_f == 2) = thetah_f;

Phi_f(S_f==1) = phil_f;
Phi_f(S_f==2) = phih_f;

U_f(1) = 0.05;
V_f(1) = U_f(1)*Theta_f(1);

for i = 2:T
  U_f(i) = U_f(i - 1)+lambda*(1 - U_f(i - 1)) - Phi_f(i - 1)*U_f(i - 1);
  V_f(i) = U_f(i)*Theta_f(i);
end

P_f_m = mean(P_f(1001:end));
P_f_sd = std(log(P_f(1001:end)));

Theta_f_m = mean(Theta_f(1001:end));
Theta_f_sd = std(log(Theta_f(1001:end)));

Phi_f_m = mean(Phi_f(1001:end));
Phi_f_sd = std(log(Phi_f(1001:end)));

V_f_m = mean(V_f(1001:end));
V_f_sd = std(log(V_f(1001:end)));

U_f_m = mean(U_f(1001:end));
U_f_sd = std(log(U_f(1001:end)));

x = linspace(1, 200, 200);
plot(x', Phi_f((T - 199):T), 'LineWidth', 2);
hold on
plot(x', U_f((T - 199):T), 'LineWidth', 2);
grid on
xlabel('Last 200 periods');
hold on

% Creating legend
legend('Job-finding rate (Blue)', 'Unemployment rate (Red)', 'Location', 'Best');

% Creating title
title({'Time series of job-finding and unemployment rates,', 'for the last 200 observations of the time period.'});

saveas(gcf, 'path\to\graphics\f_plot.png');
close(gcf);
%=======
% ANSWER
%=======
% By decreasing probability \pi from 0.96 to 0.8 (meaning 1 - \pi increases
% from 0.04 to 0.2), we see that there is an increase in the frequency of
% cycles in both the job-finding and unemployment rates. However, we see
% that the standard deviations of both time series decreases. As a result,
% we can say that the volatility of both series decreases.

% Intuitively, this is because by making the probability of switching
% between the high and low states more frequent, the two time series have
% “less time” to move to the new steady-state and achieve in a difference
% of high absolute magnitude. Instead, they will just “wiggle” between
% values whose differences are smaller. In contrast, when the probability
% 1 – \pi is smaller, that results in state changes to become more of a
% bigger “shock”, which can be characterised by the larger standard
% deviations seen in the time series of both rates.
%===========
% END ANSWER
%===========
%==========================================================================