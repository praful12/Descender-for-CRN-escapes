%traj_pts_indx = number of points in the surf plot
function [fun3d, fun3d_fft] = surf_input(arr_inp,num_spec,tot_iter,traj_pt_arr,traj_pts_indx)


fun3d = zeros(num_spec,traj_pts_indx,tot_iter);
fun3d_fft = zeros(num_spec,traj_pts_indx,tot_iter);

cum_traj = [0 ;cumsum(traj_pt_arr)];

for i = 1:tot_iter
    fun = arr_inp(cum_traj(i)+1:cum_traj(i+1),:);
    t_interp = linspace(1,traj_pt_arr(i),traj_pts_indx);
    fun = interp1([1:traj_pt_arr(i)],fun,t_interp);
    fun3d(:,:,i) = fun';
    fun3d_fft(:,:,i) = abs(fft(fun))';
end        
end    