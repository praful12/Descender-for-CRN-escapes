%Space uniform sampling

function [traj, len_curve]= space_uniform_sampling(traj_ic,traj_pts)
    diff_traj = diff(traj_ic);

    dlen_traj =(sum(diff_traj.*diff_traj,2)).^(0.5);
    dlen_traj = [0; dlen_traj];
    len_curve = cumsum(dlen_traj);
    len_curve_prime = linspace(0,len_curve(end),traj_pts);

    traj = interp1(len_curve, traj_ic,len_curve_prime);
end