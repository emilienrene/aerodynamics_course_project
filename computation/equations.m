function eq = equations()

    eq.cl2D = @(alpha) 2*pi .* deg2rad(alpha);
    
    eq.inducedAngle = @(alpha, CL) alpha - rad2deg(CL ./ (2*pi));
    
    eq.downwash = @(alpha_i, U) U .* deg2rad(alpha_i);

    eq.dw = @(dgamma, y, y0, h) -(dgamma ./ (4*pi)) .* ( ...
        1 ./ (y - y0) - ...
        (y - y0) ./ ((y - y0).^2 + (2*h).^2));

    eq.liftingLineGE = @solveLiftingLineGE;

    eq.liftingLineLLT_GE = @ solveCoupledLLT_GE;

    eq.lumpedVortex2D = @solveLumpedVortex2D;

end

function [CL, a3D, A_unit, theta] = solveLiftingLineGE(alpha_deg, AR, h, c, N)
    % Solves Prandtl's lifting-line equation with image trailing vortices
    % Inputs:
    %   alpha_deg : angle of attack vector [deg]
    %   AR        : effective aspect ratio (b / c for rectangular wing)
    %   h         : ground clearance [m]
    %   c         : chord length [m]
    %   N         : number of odd Fourier terms / collocation points (e.g. 15)

    if nargin < 5, N = 15; end

    b = AR * c;
    a_inf = 2 * pi;                         % 2D thin-airfoil lift slope
    n = 1:2:(2*N - 1);                      % Odd harmonics (symmetric wing)

    % Collocation points on half-span: theta in (0, pi/2]
    theta = (1:N)' * (pi / (2*N));          % N x 1

    % Midpoint quadrature grid for smooth image integral over [0, pi]
    M_quad = 200;
    dtheta = pi / M_quad;
    theta_q = ((1:M_quad) - 0.5) * dtheta;  % 1 x M_quad

    % Image kernel K(j, m) and cosine basis C(m, k)
    diff_cos = cos(theta_q) - cos(theta);   % N x M_quad (implicit expansion)
    K = diff_cos ./ (diff_cos.^2 + (4*h / b)^2);
    C = cos(theta_q' * n);                  % M_quad x N

    % Image integral matrix I(j, k)
    I_img = (dtheta / pi) * (K * C);        % N x N

    % Assemble linear system M * A = 1 (for unit alpha = 1 rad)
    sin_nt = sin(theta * n);                % N x N
    M_mat = (4*b / (a_inf * c)) .* sin_nt ...
          + (sin_nt ./ sin(theta) - I_img) .* n;

    % Solve for Fourier coefficients per radian of alpha
    A_unit = M_mat \ ones(N, 1);

    % 3D lift slope and lift coefficient
    a3D = pi * AR * A_unit(1);
    CL = a3D .* deg2rad(alpha_deg);
end

function CL = solveLumpedVortex2D(alpha_deg, h, c)
    % 2D flat plate in ground effect via lumped vortex (Katz & Plotkin Sec. 5.5)
    a = deg2rad(alpha_deg);
    gamma_hat = pi .* sin(a) .* ...
        (1 - (c ./ (2*h)) .* sin(a) + (c.^2 ./ (16*h.^2))) ./ ...
        (1 - (c ./ (4*h)) .* sin(a));       % Gamma / (Q_inf * c)
    CL = 2 .* gamma_hat .* (1 - (c ./ (4*pi*h)) .* gamma_hat);
end

function CL = solveCoupledLLT_GE(alpha_deg, AR, h, c, N)
    % Coupled 3D Prandtl Lifting-Line + 2D Lumped-Vortex Ground Effect (Newton solver)
    if nargin < 5, N = 15; end

    b = AR * c;
    beta = c / h;
    n = 1:2:(2*N - 1);
    theta = (1:N)' * (pi / (2*N));

    % Precompute downwash matrix W_mat (alpha_i = W_mat * A)
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

        % Initial guess from linear 3D LLT + GE
        A = M_lin \ (a_rad * ones(N, 1));

        % Newton-Raphson iterations
        for iter = 1:20
            alpha_eff = a_rad - W_mat * A;
            s = sin(alpha_eff);
            u = cos(alpha_eff);

            % 2D lumped-vortex circulation and its derivative w.r.t. alpha_eff
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

        % Integrate nonlinear sectional lift (Eq. 5.107) over half-span
        gamma_hat_h = 2 * AR * (sin_nt_h * A);
        CL(idx) = 2 * sum(gamma_hat_h .* (1 - (beta/(4*pi)) .* gamma_hat_h) ...
                          .* sin(th_h)) * dth_h;
    end
end