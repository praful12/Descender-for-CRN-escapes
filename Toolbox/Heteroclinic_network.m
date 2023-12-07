%Get the heteroclinic network and all connecting relaxation trajectories
%resampled at resamp_pt points along the relaxation trajectory.

%Input: stab_idx, sad_idx, unst_vec

num_sad = size(sad_idx,1);
resamp_pt = 4000;

%Input: pos_root_arr, stab_idx, sad_idx, unst_vec, MAK_fun, jjac_fun
%Output: hc_net, hc_traj_arr, hc_net_time, t_arr, 

hc_traj_arr = odeintegrator(@(t,q)MAK_fun_t(t,q), x_ic_arr, t_max);


hc_net = zeros(num_sad,3);
hc_net_time = zeros(num_sad,2);
x_ic_arr = zeros(num_sad,2,num_spec);
hc_traj_arr = zeros(num_sad,2,resamp_pt,num_spec);
t_arr = zeros(num_sad,2);
%get initial conditions by slightly displacing along the unstable
%eigenvalue
for i = 1:num_sad
    x_ic = pos_root_arr(sad_idx(i),:);
    [V,D] = eig(jac_fun(x_ic));
    
    for j = 1:num_spec
        if real(D(j,j))>0
            unst_vec = V(:,j)';
            break
        end
    end
    
    eps = 10^(-3);
    x_ic_arr(i,1,:) = x_ic + eps*unst_vec;
    x_ic_arr(i,2,:) = x_ic - eps*unst_vec;    
end


%Integrate MAK from each initial condition

dt = 0.5*10^(-2);
t_max = 5*10^5;
t_ran = 0:dt:t_max;
size_t = size(t_ran,2);

figure()
hold on

for i = 1:num_sad
    hc_net(i,1) = sad_idx(i);
    for j = 1:2
        "Computing time"
        tic
        conc_traj = zeros(size_t,num_spec);
        conc_traj(1,:) = x_ic_arr(i,j,:);

        for k = 2:size(t_ran,2)
            if sum(abs(MAK_fun(conc_traj(k-1,:))))>10^(-7)
                conc_traj(k,:) = conc_traj(k-1,:) + dt*MAK_fun(conc_traj(k-1,:))';
            else
                t_inter = linspace(0,t_ran(k-1),resamp_pt);
                hc_traj_arr(i,j,:,:) = interp1(t_ran(1:k-1),conc_traj(1:k-1,:),t_inter);
                break
            end
        end
        toc
        "Simulation time"
        t_ran(k)
        t_arr(i,j) = t_ran(k);
        hc_net_time(i,j) = t_ran(k);
        %use the end point to populate heteroclinic network
        end_pt = conc_traj(k-1,:);
        [a,b] = min(sum(abs(pos_root_arr(stab_idx,:)-ones(size(stab_idx,1),1)*end_pt),2));
        
        if abs(a) < 0.01
            hc_net(i,j+1) = stab_idx(b);
        end    
        %pos_root_arr(stab_idx(b),:)
        
        hc_traj_arr(i,j,1,:) = pos_root_arr(sad_idx(i),:);
        hc_traj_arr(i,j,end,:) = pos_root_arr(stab_idx(b),:);
        
        %plots
        
        plot(squeeze(hc_traj_arr(i,j,:,1)),squeeze(hc_traj_arr(i,j,:,2:end)),'Linewidth',2)
        plot(squeeze(hc_traj_arr(i,j,1,1)),squeeze(hc_traj_arr(i,j,1,2:end)),'bo','MarkerSize',10,'MarkerFaceColor','b')
        plot(squeeze(hc_traj_arr(i,j,end,1)),squeeze(hc_traj_arr(i,j,end,2:end)),'k*','MarkerSize',20)        
        
    end
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
pos_root_arr

%Also get time arrays and save
%save('..\Data\' + model_name + '_hcnet.mat', 'pos_root_arr', 'hc_net', 'hc_traj_arr', 'stab_idx', 'sad_idx','t_arr')
