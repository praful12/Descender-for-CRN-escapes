%Get the heteroclinic network and all connecting relaxation trajectories
%resampled at resamp_pt points along the relaxation trajectory.
%Input: pos_root_arr, stab_idx, sad_idx, unst_vec, MAK_fun, jjac_fun
%Output: hc_net, hc_traj_arr, hc_net_time, t_arr, 

num_sad = size(sad_idx,1);
resamp_pt = 4000;
hc_net = zeros(num_sad,3);
hc_net_time = zeros(num_sad,2);
x_ic_arr = zeros(num_sad,2,num_spec)
hc_traj_arr = cell(num_sad,2);
%t_arr = zeros(num_sad,2);
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
    x_ic_arr(i,1,:) = x_ic + eps*unst_vec
    x_ic_arr(i,2,:) = x_ic - eps*unst_vec  
end

%Integrate MAK from each initial condition using ode45
options = odeset('AbsTol',1e-10, 'RelTol', 1e-8);
dt = 0.5*10^(-2);
t_max = 5*10^5;
init_cond = zeros(1,num_spec);
t_desired = linspace(0,t_max,resamp_pt);

"Computing time"
tic
for i=1:num_sad
    hc_net(i,1) = sad_idx(i);
    for j=1:2

        init_cond = reshape(x_ic_arr(i,j,:),[],1);
        [t,sol] = ode23s(@(t,q)MAK_fun(t,q),[0 t_max], init_cond,options);

        sol_desired=interp1(t,sol,t_desired,'linear');

        hc_traj_arr{i,j} = [t_desired', sol_desired];
        
    end
end

toc
"Simulation time"

% Initialize a figure for the plot
figure;
hold on;

% Loop through each set of initial conditions and their solutions
for i = 1:num_sad
    for j=1:2
        sol = hc_traj_arr{i,j}(:, 2:end); % Extract the solution, excluding the time variable
    
        % Plot the trajectory of the solution
        plot(sol(:, 1), sol(:, 2), '-o', 'DisplayName', ['Initial Condition ' num2str(2*i+j-1)]);
    end
end

% Customize the plot (labels, title, legend, etc.)
xlabel('x-axis'); % Replace with your labels
ylabel('y-axis'); % Replace with your labels
title('Trajectories of Solutions');
legend('Location', 'best');
grid on;

% Hold off to finish plotting
hold off;

% 
% 
% ax=gca;
% ax.FontSize = 15;
% h1=legend('Relaxation','Saddle','Stable','location','best');
% %h1=legend('Stable','Saddle','location','best');
% set(h1,'Interpreter','latex');
% h1.FontSize = 15;
% 
% xlabel('$q_1$','Interpreter',"latex",'FontSize',20)
% ylabel('$q_2$','Interpreter',"latex",'FontSize',20)
% 
% hold off
% %grid on
% 
% 
% title('Proxy Heteroclinic Network','Interpreter','Latex','FontSize',20)
%saveas(gcf,'Plots/'+model_name+'_hcnet.png')




hc_net
pos_root_arr

%Also get time arrays and save
%save('..\Data\' + model_nam + '_hcnet.mat', 'pos_root_arr', 'hc_net', 'hc_traj_arr', 'stab_idx', 'sad_idx','hc_net_time')
