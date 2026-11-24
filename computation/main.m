eq = equations();

data = readtable('data/data4group20.dat', 'FileType', 'text');
data.Cl_2D = eq.cl2D(data.alpha);
data.alpha_i = eq.inducedAngle(data.alpha, data.CL);
data.w = eq.downwash(data.alpha_i, u_infinity);

p = plot_data();
p.lift(data);
p.moment(data);
