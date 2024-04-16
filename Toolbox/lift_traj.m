%Lift trajectory using numerical or analytic solver. Analytic is uncommented
%Analytic solver fails of dx = 0. For that change Eqns.
function [mom_traj, dt_traj] = lift_traj(traj, eqns_fsolve_fun,Jac_eqns_fsolve_fun,Ham_fun_qp)%,Eqns,dHamdp_fun,num_solve_end)

%get delta_x along trajectory
[traj_pts,num_spec] = size(traj);
dx_traj = diff(traj);

%Find mid-points of each segment, label it traj_prime. 
%We will assign momentum at these mid-points
traj_prime = traj(1:end-1,:)+dx_traj/2;
x_v_interp = zeros(traj_pts-1,num_spec+2);

%for i = 1:10
parfor i = 1 : traj_pts-1 % num_solve_end
    xtraj = traj_prime(i,:);
    dxtraj = dx_traj(i,:);

    P_val = eval_mom_num(xtraj,dxtraj,eqns_fsolve_fun,Jac_eqns_fsolve_fun);
    %P_val = eval_mom_anal(xtraj,dxtraj,Eqns,dHamdp_fun);

    flag = verify_const(i, P_val,xtraj,Ham_fun_qp);
    if flag == 0 
        x_v_interp(i,:) = [i P_val(2:num_spec+1) P_val(1)];
    end
end
%{
parfor i = num_solve_end + 1 : traj_pts-1
    xtraj = traj_prime(i,:);
    dxtraj = dx_traj(i,:);

    %P_val = eval_mom_num(xtraj,dxtraj,eqns_fsolve_fun,Jac_eqns_fsolve_fun);
    P_val = eval_mom_anal(xtraj,dxtraj,Eqns,dHamdp_fun);

    flag = verify_const(i, P_val,xtraj,Ham_fun_qp);
    if flag == 0 
        x_v_interp(i,:) = [i P_val(2:num_spec+1) P_val(1)];
    end
end
%}
x_v_interp( ~any(x_v_interp,2), : ) = [];
x_v_interp = [x_v_interp;
    traj_pts zeros(1,num_spec+1)];

%Display points rejected
if traj_pts  - size(x_v_interp,1) > 0
    disp("Points rejected in momentum evaluation")
    disp(traj_pts  - size(x_v_interp,1))
end
%interpolate the values over the points where the constraints are not met
%ideally there shouldn't be any - in case there are check the error by
%uncommenting conditionals in if then.

xq = 1:traj_pts;
val_traj = interp1(x_v_interp(:,1),x_v_interp(:,2:num_spec+2),xq,'linear');

%Get momentum and time assignment
mom_traj = val_traj(:,1:num_spec);
dt_traj = val_traj(:,num_spec+1);

end
