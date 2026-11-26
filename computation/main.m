%%% NACA0018 curve fitting from 3D to 2D %%%
%%% Everything runs on computation/main.m %%%

clc;
clear;
close all;

% Get all functions and variables
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
run("config\variables.m");


eq = equations();

data = readtable('data/data4group20.dat', 'FileType', 'text');
data.Cl_2D   = eq.cl2D(data.alpha);
data.alpha_i = eq.inducedAngle(data.alpha, data.CL);
data.w       = eq.downwash(data.alpha_i, u_infinity);

%% 1. Restrict fit to linear / attached regime (avoid stall > 14 deg)
linIdx = data.alpha >= -5 & data.alpha <= 10;
alpha_fit = data.alpha(linIdx);
CL_fit    = data.CL(linIdx);

%% 2. Fit (AR, h) using nonlinear least squares (fminsearch)
N_modes = 20;
function err = lossCoupled(p, alpha_fit, CL_fit, c, N)
    AR = p(1); h = p(2);
    if AR < 0.5 || AR > 20 || h < 0.1*c || h > 5*c
        err = 1e9;
        return;
    end
    CL_pred = eq.liftingLineGE(alpha_fit, AR, h, c, N);
    err = sum((CL_pred - CL_fit).^2);
end
%lossFn = @(p) lossCoupled(alpha_fit, p(1), p(2), c, N_modes) ;
lossFn = @(p) sum((eq.liftingLineGE(alpha_fit, p(1), p(2), c, N_modes) - CL_fit).^2);

% Initial guess: [AR_guess, h_guess (m)]
p0 = [4.0, 0.15];
opts = optimset('Display', 'off');
p_opt = fminsearch(lossFn, p0, opts);

AR_opt = p_opt(1);
h_opt  = p_opt(2);
fprintf('Fitted effective AR = %.3f, h = %.4f m (h/c = %.2f)\n', ...
    AR_opt, h_opt, h_opt / c);

%% 3. Evaluate models and plot comparison
data.CL_LLT_GE = eq.liftingLineGE(data.alpha, AR_opt, h_opt, c, N_modes);
data.CL_2D_GE  = eq.lumpedVortex2D(data.alpha, h_opt, c);

p = plot_data();
p.lift(data);
p.moment(data);

% Quick overlay plot of experimental data vs. lifting-line ground effect
figure; hold on; grid on;
errorbar(data.alpha, data.CL, data.CLstd, 'ko', 'DisplayName', 'Experiment');
plot(data.alpha, data.Cl_2D, '--', 'DisplayName', '2D Thin Airfoil');
plot(data.alpha, data.CL_LLT_GE, 'r-', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('3D LLT + GE (AR=%.2f, h=%.3fm)', AR_opt, h_opt));
plot(data.alpha, data.CL_2D_GE, 'b-.', 'LineWidth', 1.2, ...
    'DisplayName', '2D Lumped Vortex + GE');
xlabel('\alpha [deg]'); ylabel('C_L');
xlim([-6 16]); ylim([-0.5 1.5]);
legend('Location', 'northwest');