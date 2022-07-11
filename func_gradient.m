%The momentum assignment is done at the middle of each segment.
%However in PS_traj, the momentum assignment is at the same point as the
%position, thus we take the mean of the neighboring momentum to define
%mom_prime which then is used for PS_traj.


function [PS_traj, delta_x, t_traj] = func_gradient(traj, mom_traj, dt_traj, dHamdq_fun,num_smooth)

mom_traj_ic = mom_traj;
dt_traj_ic = dt_traj;
[traj_pts,num_spec] = size(traj);
for j = 1:num_smooth
    dt_traj = smooth(dt_traj);
    for i = 1:num_spec
        mom_traj(:,i) = smooth(mom_traj(:,i));
    end
end


%{
figure()
plot(mom_traj,'--')
hold on
plot(mom_traj_ic)
hold off

figure()
plot(dt_traj,'--')
hold on
plot(dt_traj_ic)
hold off
%}

%Get momentum at the end points. Recall that the momentum assignment was 
%done at the middle of the segments.
mom_prime = zeros(size(traj));
diff_mom = diff(mom_traj);
mom_prime(2:end-1,:) = mom_traj(1:end-2,:) + 1/2*diff_mom(1:end-1,:);
%diff_mom = diff(mom_prime);

%Time taken to go from mid point i to i+1 is (dt(i)+dt(i+1))/2
dt_prime = 1/2*(dt_traj(1:end-2)+dt_traj(2:end-1));

%delta_x and p_dot
p_dot = zeros(size(traj));
delta_x = zeros(size(traj));

%Output
PS_traj = horzcat(traj,mom_prime);
t_traj =  [0;cumsum(dt_traj(1:end-1))];    
p_dot(2:end-1,:) = diff_mom(1:end-1,:)./dt_prime;
for i = 2:traj_pts-1
    delta_x(i,:) = p_dot(i,:) + dHamdq_fun(PS_traj(i,:));
end

%{
figure()
plot(delta_x)
title('delta-x')

figure()
plot(p_dot)
title('p-dot')

figure()
plot(dHamdq_fun(PS_traj(:,:)))
title('dHdq')
%}

end