%A function to pick the step size

function [eps, action_S_new, PS_traj, t_traj, delta_x, S_traj] = eps_picker(traj_or,eqns_fsolve_fun,Jac_eqns_fsolve_fun,dHamdq_fun,Ham_fun_qp,delta_x_or,action_S,eps_ic,num_smooth,eps_min)

eps = eps_ic;
delta_S = 1;
eps = eps*2;
%Calculate action, if delta_S > 0, then half the epsilon.
while delta_S > 0  && eps > eps_min %10^(-4)
    eps = eps*0.5;
    traj = traj_or + eps*delta_x_or;
    [PS_traj, t_traj, delta_x] = traj_2_PS_traj_func(traj,eqns_fsolve_fun,Jac_eqns_fsolve_fun,dHamdq_fun,Ham_fun_qp,num_smooth);
    [action_S_new, S_traj] = action_PS_traj(PS_traj);
    delta_S = action_S_new - action_S
    %eps_ct = eps_ct + 1;    
end

if  eps <= eps_min
    eps = 0;
    [PS_traj, t_traj, delta_x] = traj_2_PS_traj_func(traj_or,eqns_fsolve_fun,...
        Jac_eqns_fsolve_fun,dHamdq_fun,Ham_fun_qp,num_smooth);
    [action_S_new, S_traj] = action_PS_traj(PS_traj); 
    action_S_new = action_S;
end

end