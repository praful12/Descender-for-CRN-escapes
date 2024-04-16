function traj_ic = straight_line_IC(traj_pt,pos_root_arr,start_idx,end_idx)

p_o_int = pos_root_arr(start_idx,:);
p_o_int = [p_o_int; pos_root_arr(end_idx,:)];

t_fin = 1;
t_i = 0 : t_fin/(traj_pt-1) : t_fin;
traj_ic = t_i'*(p_o_int(2,:))+(1-t_i')*(p_o_int(1,:));
