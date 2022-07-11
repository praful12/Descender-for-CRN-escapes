%Function that will pick out a trajectory from the array saved in
%update_descender
function traj_idx = picker_arr_idx(arr, traj_pt_arr, idx)
    
    cum_traj = [0 ;cumsum(traj_pt_arr)];
    traj_idx = arr(cum_traj(idx)+1:cum_traj(idx+1),:);
end