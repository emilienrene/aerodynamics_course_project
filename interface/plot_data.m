function p = plot_data()
p.lift = @plotLift;
p.moment = @plotMoment;
end

function plotLift(data)
figure;
plot(data.alpha, data.CL, 'o-', ...
    data.alpha, data.Cl_2D, 'o-');
xlabel('\alpha (°)');
ylabel('C_L');
title('Lift coefficient versus angle of attack');
legend('3D', '2D', 'Location', 'best');
grid on;
end

function plotMoment(data)
figure;
plot(data.alpha, data.CM_LE, 'o-');
xlabel('\alpha (°)');
ylabel('C_{M,LE}');
title('Moment coefficient versus angle of attack');
grid on;
end