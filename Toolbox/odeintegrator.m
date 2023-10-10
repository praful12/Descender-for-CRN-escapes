%Integrator

function [solutions] = odeintegrator(difeq, init_cond)

solutions = cell(size(init_cond,1),1);

%Integrate MAK from each initial condition using ode23s
options = odeset('AbsTol',1e-10, 'RelTol', 1e-8);
t_max = 5*10^5;

%Number of data points to resample time
resamp_pt=4000;
t_desired = linspace(0,t_max,resamp_pt);

"Computing time"
tic
parfor i=1:size(init_cond,1)
    y0= init_cond(i,:);

    [t,sol] = ode23s(difeq,[0 t_max],y0,options);

    sol_desired=interp1(t,sol,t_desired,'linear');

    solutions{i} = [t_desired', sol_desired];
end

toc
"Simulation time"