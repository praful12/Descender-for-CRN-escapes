
%Forward Euler integration of Hamilton's equations of motion.
%Give input as target point, dt and time of integration
%along with starting point and gradients. The sign of dt and t_max will
%determine forward or backward integration.
%Note: We integrate Hamiltonian trajectory from some initial conditions to as far as
%distance to target point starts increasing rather than decreasing

function [PS_traj, t_traj, dist_curr] = hamEoM_f_Euler(x_ic,x_targ,dt,t_max,dHamdq_fun,dHamdp_fun)

num_spec = size(x_ic,2)/2;
t_ran = 0:dt:t_max;
int_HamNP = zeros(size(t_ran,2),num_spec*2);
int_HamNP(1,:) = x_ic;

i = 1;
dist_curr = sum((x_targ-int_HamNP(i,1:num_spec)).^2,'all');
dist_old = dist_curr;
min_dist = 1000;

while dist_old - dist_curr > - 0.000001 && i < size(t_ran,2)
    i = i+1;
    dist_old = dist_curr;
    
    %Integrate HamEoM
    int_HamNP(i,1:num_spec) = int_HamNP(i-1,1:num_spec) + dt*dHamdp_fun(int_HamNP(i-1,:));
    int_HamNP(i,num_spec+1:2*num_spec) = int_HamNP(i-1,num_spec+1:2*num_spec) - dt*dHamdq_fun(int_HamNP(i-1,:));    
    
    %Update distance
    dist_curr = sum((x_targ-int_HamNP(i,1:num_spec)).^2,'all');
    if min_dist > dist_curr
        min_dist = dist_curr;
    end
end

PS_traj = int_HamNP(1:i,:);
t_traj = t_ran(1,1:i);
dist_curr = min_dist^0.5;
end