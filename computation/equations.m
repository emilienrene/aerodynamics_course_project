function eq = equations()

    eq.cl2D = @(alpha) 2*pi .* deg2rad(alpha);
    
    eq.inducedAngle = @(alpha, CL) alpha - rad2deg(CL ./ (2*pi));
    
    eq.downwash = @(alpha_i, U) U .* deg2rad(alpha_i);

    eq.gamma = [];

    eq.dw = @(dgamma, y, y0, h) dgamma .* ( ...
        -1 ./ (4*pi .* (y - y0)) + ...
        (y - y0) ./ (2*h.^2 + (y - y0).^2));

end