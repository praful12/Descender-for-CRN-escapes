
function peak_freq = wigner_dist_analysis_noplots(delta_x)
[traj_pts,num_spec] = size(delta_x);
pos_plots = zeros(num_spec,2*traj_pts);
freq_plots = zeros(num_spec,traj_pts);
peak_freq_idx = zeros(num_spec,1);

for i = 1:num_spec
    [d,f,t] = wvd(delta_x(:,i));
    pos_plots(i,:) = sum(d,1);
    freq_plots(i,:) = sum(d,2);
    [w,p] = findpeaks(freq_plots(i,:));
    peak_freq_idx(i) = p(1);
end

peak_freq = min(peak_freq_idx)+1;

end