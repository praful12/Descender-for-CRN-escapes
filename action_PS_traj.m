%------action along trajectory

function [action_S_new, S_traj] = action_PS_traj(PS_traj)
    num_spec = size(PS_traj,2)/2;
    traj_pts = size(PS_traj,1);
    
    traj = PS_traj(:,1:num_spec);
    mom_traj = PS_traj(:,1+num_spec:2*num_spec);
    S_traj = zeros(traj_pts,1);
    dS = 0;
    S_0 = 0;
    dq = diff(traj);
    dq = [dq ; zeros(1,num_spec)];

    for i = 2:traj_pts
        dS = mom_traj(i,:)*dq(i,:).';
        S_traj(i) = S_traj(i-1) + dS ;
    end
    action_S_new = S_traj(i);
end
