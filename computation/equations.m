function eq = equations()

    eq.cl2D = @(alpha) 2*pi .* deg2rad(alpha);
    
    eq.inducedAngle = @(alpha, CL) alpha - rad2deg(CL ./ (2*pi));
    
    eq.downwash = @(alpha_i, U) U .* deg2rad(alpha_i);

    eq.dw = @(dgamma, y, y0, h) -(dgamma ./ (4*pi)) .* ( ...
        1 ./ (y - y0) - ...
        (y - y0) ./ ((y - y0).^2 + (2*h).^2));

    % Coupled 3D Prandtl Lifting-Line + 2D Lumped-Vortex Ground Effect
    eq.coupledLLT_GE = @solveCoupledLLT_GE;

    % Center of pressure x_cp / c = -CM_LE / CL (masked where |CL| is near 0)
    eq.xcp = @(CL, CM_LE) -CM_LE ./ CL;

end

function [CL, dist] = solveCoupledLLT_GE(alpha_deg, AR, h, c, U_inf, N)
    if nargin < 5 || isempty(U_inf), U_inf = 1.0; end
    if nargin < 6 || isempty(N),     N = 15;      end

    b = AR * c;
    beta = c / h;
    n = 1:2:(2*N - 1);                      % Odd Fourier modes
    theta = (1:N)' * (pi / (2*N));          % Half-span collocation points (0, pi/2]

    % Downwash matrix W_mat such that alpha_i = W_mat * A
    M_quad = 200;
    dtheta = pi / M_quad;
    theta_q = ((1:M_quad) - 0.5) * dtheta;

    diff_cos = cos(theta_q) - cos(theta);
    K = diff_cos ./ (diff_cos.^2 + (4*h / b)^2);
    I_img = (dtheta / pi) * (K * cos(theta_q' * n));

    sin_nt = sin(theta * n);
    W_mat  = (sin_nt ./ sin(theta) - I_img) .* n;
    M_lin  = (2*AR / pi) .* sin_nt + W_mat;

    % Half-span quadrature grid [0, pi/2] for CL integration
    dth_h = (pi/2) / 100;
    th_h  = ((1:100) - 0.5)' * dth_h;
    sin_nt_h = sin(th_h * n);

    CL = zeros(size(alpha_deg));

    for idx = 1:numel(alpha_deg)
        a_rad = deg2rad(alpha_deg(idx));

        % Start from linear 3D LLT + GE solution
        A = M_lin \ (a_rad * ones(N, 1));

        % Newton-Raphson coupling loop
        for iter = 1:20
            alpha_eff = a_rad - W_mat * A;
            s = sin(alpha_eff);
            u = cos(alpha_eff);

            % 2D lumped-vortex circulation and derivative w.r.t. alpha_eff
            num  = (1 + beta^2/16) .* s - (beta/2) .* s.^2;
            den  = 1 - (beta/4) .* s;
            dnum = (1 + beta^2/16) - beta .* s;
            dden = -beta/4;

            gamma_2D  = pi .* (num ./ den);
            dgamma_da = pi .* ((dnum .* den - num .* dden) ./ (den.^2)) .* u;

            % Residual and Jacobian
            R = 2*AR * (sin_nt * A) - gamma_2D;
            J = 2*AR * sin_nt + (dgamma_da .* W_mat);

            dA = -J \ R;
            A  = A + dA;

            if norm(dA) < 1e-10
                break;
            end
        end

        % Integrate nonlinear sectional lift over half-span
        gamma_hat_h = 2 * AR * (sin_nt_h * A);
        CL(idx) = 2 * sum(gamma_hat_h .* (1 - (beta/(4*pi)) .* gamma_hat_h) ...
                          .* sin(th_h)) * dth_h;
    end

    % Build full-span distributions for the last evaluated angle of attack
    if nargout > 1
        y_half = (b/2) * cos(theta);                        % from +b/2 toward 0
        Gamma_half = 2 * b * U_inf * (sin_nt * A);          % [m^2/s]
        alpha_eff_deg = rad2deg(a_rad - W_mat * A);         % [deg]

        % Mirror half-span to full span [-b/2, +b/2] and include wingtips (Gamma = 0)
        dist.y = [-b/2; -y_half; flipud(y_half(1:end-1)); b/2];
        dist.Gamma = [0; Gamma_half; flipud(Gamma_half(1:end-1)); 0];
        dist.alpha_eff = [0; alpha_eff_deg; flipud(alpha_eff_deg(1:end-1)); 0];
    end
end