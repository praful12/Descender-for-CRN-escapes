%Optimal resampling
%If new coordinates are s, then ds = dx*1/(|delta_x|+0.1)
%So new ds is inversely proportional to magnitude of delta_x

function [traj_prime,len_curve_prime] = optimal_resampling(traj,delta_x,f0)

[traj_pts,num_spec] = size(traj);

diff_traj = diff(traj);
dlen_traj =(sum(diff_traj.*diff_traj,2)).^(0.5);
dlen_traj = [0; dlen_traj];
len_curve = cumsum(dlen_traj);
%this should be a straight line
%figure()
%plot(len_curve)

%Number of segments for resampling and indices of each segment
%For this code, there can only be 4 segments. 
seg_4_resamp = 4;               
idx_seg = ones(seg_4_resamp,2);
for i = 1:seg_4_resamp-1
    end_idx = find(len_curve - i/seg_4_resamp*len_curve(end)>0,1);
    idx_seg(i,2) = end_idx-1;    idx_seg(i+1,1) = end_idx;
end
idx_seg(end,2) = size(len_curve,1);

%Get number of points in eachs segment : b_arr
b_arr = diff(idx_seg,1,2)+1;
if sum(b_arr) - traj_pts ~= 0
    "Error in indexing"
end

%{
%plot the different segments in the trajectory
figure()
hold on
for i = 1:seg_4_resamp
    plot(traj(idx_seg(i,1):idx_seg(i,2),1),traj(idx_seg(i,1):idx_seg(i,2),2),'DisplayName',string(i))
end
%xlim([2 10])
legend 
hold off

%plot the different delta_x segments in the trajectory
figure()
hold on
for i = 1:seg_4_resamp
    plot([idx_seg(i,1):idx_seg(i,2)],delta_x(idx_seg(i,1):idx_seg(i,2),:),'DisplayName',string(i))
end
%xlim([2 10])
%legend 
hold off
%}
%Get delta_x, their fft, and reverse cumsum of |fft|^2 for the segments
ft_delta_x_arr = cell(seg_4_resamp,1);
for i = 1:seg_4_resamp
    del_x = delta_x(idx_seg(i,1):idx_seg(i,2),:);
    ft_delx = fft(del_x)/b_arr(i);
    ft_delx = abs(ft_delx).^2;
    ft_delx = ft_delx(1:floor(size(ft_delx,1)/2),:);
    cum_ftdx = cumsum(ft_delx,"reverse");
    ft_delta_x_arr{i} =cum_ftdx;
end

%ft_delta_x_arr

%{
%Plot the cumsum to predict how the result of sampling should look like
figure()
hold on
for i = 1:seg_4_resamp
    cum_arr = ft_delta_x_arr{i};
    %size(cum_arr)
    plot(cum_arr(:,2),'DisplayName',string(i))
%    plot(cum_arr(:,2),'DisplayName',string(i))
end
xlim([2 10])
ylabel('C(f)')
xlabel('f')
legend 
hold off
%}
fun_2min = @(c) energy_4_resamp(c,b_arr,f0+1,ft_delta_x_arr);

%poor man's optimization
min_each = 300;
step_sz = 50;
max_each = (traj_pts-min_each)/3 - mod((traj_pts-min_each)/3,step_sz);
samp_vec = [min_each:step_sz:max_each];
samp_steps = size(samp_vec,2);

func_val = zeros((samp_steps)^3,1);
vec_val = zeros((samp_steps)^3,4);

%Change the number of loops if number of segments are changed. 
%#loops = seg_4_resamp-1
ct =  0;
for i = 1:samp_steps
    for j = 1:samp_steps
        for k = 1:samp_steps
            ct = ct+1;
            v1 = samp_vec(i); v2 = samp_vec(j); v3 = samp_vec(k); 
            vec_val(ct,:) = [traj_pts-v1-v2-v3 v1 v2 v3];
            func_val(ct) = fun_2min(vec_val(ct,:));
        end
    end
end

[a,b] = min(func_val);
c_opt = vec_val(b,:)
%figure()
%plot(func_val)

%Now redistribute points along the trajectory
len_curve_prime = zeros(traj_pts,1);
idx_curr = 1;
for i = 1:seg_4_resamp-1
    idx_end =  idx_curr+c_opt(i)-1; 
    samp_bw = linspace(len_curve(idx_seg(i,1)),len_curve(idx_seg(i,2)+1),c_opt(i)+1);
    len_curve_prime(idx_curr:idx_end,1) = samp_bw(1:c_opt(i));
    idx_curr = idx_end+1;
end
samp_bw = linspace(len_curve(idx_seg(end,1)),len_curve(end),c_opt(end));
len_curve_prime(idx_curr:end,1) = samp_bw(1:c_opt(end));

%Smoothen the effect of discrete jumps
len_curve_diff = len_curve_prime - linspace(len_curve_prime(1),len_curve_prime(end),traj_pts)';
len_curve_diff_s =  smooth_traj_butter(len_curve_diff,4);
len_curve_prime = len_curve_diff_s + linspace(len_curve_prime(1),len_curve_prime(end),traj_pts)';

%len_curve_prime = smooth(smooth(len_curve_prime));%,'sgolay',4);
%{
figure()
plot(diff(len_curve_prime))
hold on
plot(diff(len_curve),'--k')
%plot(diff(len_curve))
hold off
%}

traj_prime = interp1(len_curve,traj,len_curve_prime);

end
