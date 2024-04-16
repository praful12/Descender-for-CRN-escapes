%This function estimates the index on a PS_traj, along which when
%integrated Hamilton's equations of motion (HamEoM) approach closest to the
%end points.

function a = ham_closest_approach_noplot(PS_traj,S_traj,t_traj,...
    dHamdp_fun, dHamdq_fun, save_plot_name,save_plot_flag,plotnam)

%Initialize
[traj_pts,num_spec] = size(PS_traj);
num_spec = num_spec/2;
traj = PS_traj(:,1:num_spec);
x_targ_start = traj(1,:);
x_targ_end = traj(end,:);
dt = min(diff(t_traj));

%Coarse search -  every 100 points
pt_arr = [1:100:traj_pts];
dist_min_arr = zeros(size(pt_arr,2),2);

for i = 1:size(pt_arr,2)
    pt_idx = pt_arr(i);
    t_for = 2*(t_traj(end)-t_traj(pt_idx));
    t_back = 2*(t_traj(1)-t_traj(pt_idx));
    x_ic = PS_traj(pt_idx,:);
    [Ham_traj_for, ~, dist_min_arr(i,1)] = hamEoM_f_Euler(x_ic,x_targ_end,dt,t_for,dHamdq_fun,dHamdp_fun);
    [Ham_traj_back, ~, dist_min_arr(i,2)] = hamEoM_f_Euler(x_ic,x_targ_start,-dt,t_back,dHamdq_fun,dHamdp_fun);
end
[a,b] = min(sum(dist_min_arr,2));

%Finer search - every 10 points

pt_arr = [pt_arr(max(1,b-1)):10:pt_arr(min(size(pt_arr,2),b+1))];
dist_min_arr = zeros(size(pt_arr,2),2);

for i = 1:size(pt_arr,2)
    pt_idx = pt_arr(i);
    t_for = 2*(t_traj(end)-t_traj(pt_idx));
    t_back = 2*(t_traj(1)-t_traj(pt_idx));
    x_ic = PS_traj(pt_idx,:);
    [Ham_traj_for, ~, dist_min_arr(i,1)] = hamEoM_f_Euler(x_ic,x_targ_end,dt,t_for,dHamdq_fun,dHamdp_fun);
    [Ham_traj_back, ~, dist_min_arr(i,2)] = hamEoM_f_Euler(x_ic,x_targ_start,-dt,t_back,dHamdq_fun,dHamdp_fun);
end
[a,b] = min(sum(dist_min_arr,2));


end