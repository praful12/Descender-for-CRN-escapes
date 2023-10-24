%Get the heteroclinic network and all connecting relaxation trajectories
%resampled at resamp_pt points along the relaxation trajectory.
%Input: pos_root_arr, stab_idx, sad_idx, unst_vec, MAK_fun, jjac_fun
%Output: hc_net, hc_traj_arr, hc_net_time, t_arr, 

hc_traj_arr = odeintegrator(@(t,q)MAK_fun(t,q), x_ic_arr, t_max);

num_sad = size(sad_idx,1);
hc_net = zeros(num_sad,3);

% Initialize a figure for the plot
figure;
hold on;

% Loop through each set of initial conditions and their solutions
for i = 1:size(x_ic_arr,1)
    saddle_pt = floor((i-1)/2)+1;
    direction = mod(i-1,2)+1;

    hc_net(saddle_pt,1) = sad_idx(saddle_pt);

    sol = hc_traj_arr{i}(:, 2:end); % Extract the solution, excluding the time variable

    end_pt = hc_traj_arr{i}(end,2:end);
    [a,b] = min(sum(abs(pos_root_arr(stab_idx,:)-ones(size(stab_idx,1),1)*end_pt),2));
    
    if abs(a) < 0.01
        hc_net(saddle_pt,direction+1) = stab_idx(b);
    end    
    %pos_root_arr(stab_idx(b),:)
    
    hc_traj_arr{i}(1,2:end) = pos_root_arr(sad_idx(saddle_pt),:);
    hc_traj_arr{i}(end,2:end) = pos_root_arr(stab_idx(b),:);
    
    %plots
    
    plot(sol(:,1),sol(:,2:end),'Linewidth',2)
    plot(hc_traj_arr{i}(1,2),hc_traj_arr{i}(1,3:end),'bo','MarkerSize',10,'MarkerFaceColor','b')
    plot(hc_traj_arr{i}(end,2),hc_traj_arr{i}(end,3:end),'k*','MarkerSize',20)        
end

ax=gca;
ax.FontSize = 15;
h1=legend('Relaxation','Saddle','Stable','location','best');
%h1=legend('Stable','Saddle','location','best');
set(h1,'Interpreter','latex');
h1.FontSize = 15;

xlabel('$q_1$','Interpreter',"latex",'FontSize',20)
ylabel('$q_2$','Interpreter',"latex",'FontSize',20)

hold off
%grid on


title('Proxy Heteroclinic Network','Interpreter','Latex','FontSize',20)
%saveas(gcf,'Plots/'+model_name+'_hcnet.png')

hc_net
