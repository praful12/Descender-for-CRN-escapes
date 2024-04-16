%number of smoothenings
num_smooth = 45;

%Initialize descender
[PS_traj, t_traj, delta_x] = traj_2_PS_traj_func(traj...
    ,eqns_fsolve_fun,Jac_eqns_fsolve_fun,dHamdq_fun,Ham_fun_qp,num_smooth);

[action_S, S_traj] = action_PS_traj(PS_traj);
f0 = wigner_dist_analysis_noplots(delta_x);
delta_x_s =  smooth_traj_butter(delta_x,f0); 


%Save to file
traj_pt_arr = [traj_pts];
PS_arr = [PS_traj];
deltax_arr = [delta_x];
deltax_lp_arr = [delta_x_s];
time_arr = [t_traj];
S_traj_arr = [S_traj];
S_arr = [action_S]
f0_arr = [f0];
eps_arr = [];
a_arr = [];

%s_resamp_arr = [len_curve_prime];