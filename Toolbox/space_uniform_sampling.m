%Space uniform sampling

function [traj, len_curve]= space_uniform_sampling(traj_ic,traj_pts)
    diff_traj = diff(traj_ic);
    
    %Sometimes the initial array doesnt budge. So we have to remove the
    %extra points in order to linearly interpolate
    dlen_traj =(sum(diff_traj.*diff_traj,2)).^(0.5);
    I = find(dlen_traj~=0);
    traj_ic = [traj_ic(1,:); traj_ic(I+1,:)];

    diff_traj = diff(traj_ic);    
    dlen_traj =(sum(diff_traj.*diff_traj,2)).^(0.5);
    dlen_traj = [0; dlen_traj];
    len_curve = cumsum(dlen_traj);
    len_curve_prime = linspace(0,len_curve(end),traj_pts);
    
    traj = interp1(len_curve, traj_ic,len_curve_prime);
end