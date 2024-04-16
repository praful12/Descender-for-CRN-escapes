function [PS_traj, t_traj, delta_x] = traj_2_PS_traj_func(traj,eqns_fsolve_fun,Jac_eqns_fsolve_fun,dHamdq_fun,Ham_fun_qp,num_smooth)

%lift trajectory
[mom_traj, dt_traj] = lift_traj(traj, eqns_fsolve_fun,Jac_eqns_fsolve_fun,Ham_fun_qp);


%get PS_traj, gradient and t_traj
[PS_traj, delta_x, t_traj] = func_gradient(traj, mom_traj, dt_traj, dHamdq_fun,num_smooth);


end
