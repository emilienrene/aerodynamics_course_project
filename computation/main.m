%%% NACA0018 curve fitting from 3D to 2D %%%
%%% Everything runs on computation/main.m %%%

clc;
clear;
close all;

% Get all functions and variables
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
run("config\variables.m");

if ~exist('c', 'var'), c = 0.15; end
if ~exist('u_infinity', 'var'), u_infinity = 1.0; end

eq = equations();

data = readtable('data/data4group20.dat', 'FileType', 'text');
data.Cl_2D   = eq.cl2D(data.alpha);
data.alpha_i = eq.inducedAngle(data.alpha, data.CL);
data.w       = eq.downwash(data.alpha_i, u_infinity);

%% 1. Fit (AR, h) in the attached-flow regime (-5 to 10 deg)
linIdx    = data.alpha >= -5 & data.alpha <= 10;
alpha_fit = data.alpha(linIdx);
CL_fit    = data.CL(linIdx);

N_modes = 15;
lossFn  = @(p) lossCoupled(p, eq, alpha_fit, CL_fit, c, u_infinity, N_modes);

p0    = [4.0, 0.10];    % Initial guess: [AR, h (m)]
opts  = optimset('Display', 'off');
p_opt = fminsearch(lossFn, p0, opts);

AR_opt = p_opt(1);
h_opt  = p_opt(2);
fprintf('Fitted effective AR = %.3f, h = %.4f m (h/c = %.2f)\n', ...
    AR_opt, h_opt, h_opt / c);

%% 2. Evaluate fitted curve and spanwise distributions at a chosen alpha
data.CL_coupled = eq.coupledLLT_GE(data.alpha, AR_opt, h_opt, c, u_infinity, N_modes);

alpha_plot = 8; % Angle of attack [deg] for spanwise distribution plots
[~, dist] = eq.coupledLLT_GE(alpha_plot, AR_opt, h_opt, c, u_infinity, N_modes);

%% 3. Plots
p = plot_data();
p.lift(data);
p.moment(data);

% Plot 1: Experimental CL vs Coupled 3D LLT + GE fit
figure; hold on; grid on;
errorbar(data.alpha, data.CL, data.CLstd, 'ko', 'DisplayName', 'Experiment');
plot(data.alpha, data.Cl_2D, '--', 'DisplayName', '2D Thin Airfoil');
plot(data.alpha, data.CL_coupled, 'r-', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('Coupled LLT + GE (AR=%.2f, h/c=%.2f)', AR_opt, h_opt/c));
xlabel('\alpha [deg]'); ylabel('C_L');
xlim([-6 16]); ylim([-0.5 1.5]);
legend('Location', 'northwest');

% Plot 2: Spanwise alpha_eff(y) and Gamma(y) distributions
figure;
subplot(2, 1, 1); hold on; grid on;
plot(dist.y, dist.alpha_eff, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
yline(alpha_plot, '--k', sprintf('Geometric \\alpha = %g^\\circ', alpha_plot));
xlabel('Spanwise position y [m]');
ylabel('\alpha_{eff} [deg]');
title(sprintf('Spanwise Distributions at \\alpha = %g^\\circ (AR = %.2f, h/c = %.2f)', ...
    alpha_plot, AR_opt, h_opt/c));

subplot(2, 1, 2); hold on; grid on;
plot(dist.y, dist.Gamma, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Spanwise position y [m]');
ylabel('\Gamma(y) [m^2/s]');

%% Local helper function for bounded least squares
function err = lossCoupled(p, eq, alpha_fit, CL_fit, c, U_inf, N)
    AR = p(1);
    h  = p(2);
    if AR < 0.5 || AR > 20 || h < 0.1*c || h > 5*c
        err = 1e9;
        return;
    end
    CL_pred = eq.coupledLLT_GE(alpha_fit, AR, h, c, U_inf, N);
    err = sum((CL_pred - CL_fit).^2);
end