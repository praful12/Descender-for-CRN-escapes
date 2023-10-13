
%----- update and save iteration

action_S = action_S_new;

traj_pt_arr = [traj_pt_arr ; traj_pts];
PS_arr = [PS_arr; PS_traj];
deltax_arr = [deltax_arr; delta_x];
deltax_lp_arr = [deltax_lp_arr; delta_x_s];
time_arr = [time_arr; t_traj];
S_traj_arr = [S_traj_arr; S_traj];
S_arr = [S_arr; action_S];
f0_arr = [f0_arr; f0];
eps_arr = [eps_arr; eps_used];
s_resamp_arr = [s_resamp_arr;len_curve_prime];
a_arr = [a_arr; a_min];