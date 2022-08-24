%This function estimates the index on a PS_traj, along which when
%integrated Hamilton's equations of motion (HamEoM) approach closest to the
%end points.

function [a,Ham_traj_for, Ham_traj_back] = Ham_closest_approach(PS_traj,S_traj,t_traj,...
    dHamdp_fun, dHamdq_fun, save_plot_name,save_plot_flag,plotnam)

%Initialize
[traj_pts,num_spec] = size(PS_traj);
num_spec = num_spec/2;
traj = PS_traj(:,1:num_spec);
x_targ_start = traj(1,:);
x_targ_end = traj(end,:);
dt = min(diff(t_traj));

%Coarse search -  every 100 points
pt_arr = [1:50:traj_pts];
dist_min_arr = zeros(size(pt_arr,2),2);

for i = 1:size(pt_arr,2)
    pt_idx = pt_arr(i);
    t_for = 2*(t_traj(end)-t_traj(pt_idx));
    t_back = 2*(t_traj(1)-t_traj(pt_idx));
    x_ic = PS_traj(pt_idx,:);
    [Ham_traj_for, ~, dist_min_arr(i,1)] = HamEoM_f_Euler(x_ic,x_targ_end,dt,t_for,dHamdq_fun,dHamdp_fun);
    [Ham_traj_back, ~, dist_min_arr(i,2)] = HamEoM_f_Euler(x_ic,x_targ_start,-dt,t_back,dHamdq_fun,dHamdp_fun);
end
[a,b] = min(sum(dist_min_arr,2));

%Finer search - every 10 points
pt_arr = [pt_arr(max(1,b-1)):10:pt_arr(min(size(pt_arr,2),b+1))];
dist_min_arr = zeros(size(pt_arr,2),2);

for i = 1:size(pt_arr,2)
    pt_idx = pt_arr(i);
    t_for = 2*(t_traj(end)-t_traj(pt_idx));
    t_back = 2*(t_traj(1)-t_traj(pt_idx));
    x_ic = PS_traj(pt_idx,:);
    [Ham_traj_for, ~, dist_min_arr(i,1)] = HamEoM_f_Euler(x_ic,x_targ_end,dt,t_for,dHamdq_fun,dHamdp_fun);
    [Ham_traj_back, ~, dist_min_arr(i,2)] = HamEoM_f_Euler(x_ic,x_targ_start,-dt,t_back,dHamdq_fun,dHamdp_fun);
end
[a,b] = min(sum(dist_min_arr,2));

%Final trajectory
pt_idx = pt_arr(b);
t_for = 2*(t_traj(end)-t_traj(pt_idx));
t_back = 2*(t_traj(1)-t_traj(pt_idx));
x_ic = PS_traj(pt_idx,:);
[Ham_traj_for, t_f, dist_curr_for] = HamEoM_f_Euler(x_ic,x_targ_end,dt,t_for,dHamdq_fun,dHamdp_fun);
t_f = t_f + t_traj(pt_idx);
[Ham_traj_back, t_b, dist_curr_back] = HamEoM_f_Euler(x_ic,x_targ_start,-dt,t_back,dHamdq_fun,dHamdp_fun);
t_b = t_b + t_traj(pt_idx);


%Plot
figure()
subplot(1,2,1)
plot(traj(:,1),traj(:,2),'c-.','LineWidth',2)
hold on
plot(Ham_traj_for(:,1),Ham_traj_for(:,2),'LineWidth',2)
plot(Ham_traj_back(:,1),Ham_traj_back(:,2),'b','LineWidth',2)
plot(traj(:,1),traj(:,3:end),'c-.','LineWidth',2)
plot(Ham_traj_for(:,1),Ham_traj_for(:,3:num_spec),'LineWidth',2)
plot(Ham_traj_back(:,1),Ham_traj_back(:,3:num_spec),'b','LineWidth',2)
plot(x_ic(:,1),x_ic(2:num_spec),'gs','MarkerSize',15)
plot(x_targ_start(:,1),x_targ_start(:,2:end),'k*','MarkerSize',15)
plot(x_targ_end(:,1),x_targ_end(:,2:end),'bo','MarkerSize',10,'MarkerFaceColor','b')
%grid on
hold off

%legend('Descender','Fwd HamEoM','Bwd HamEoM','location','best')
xlabel('$q_1$','Interpreter','Latex','Fontsize',20)
ylabel('$q_2$ onwards','Interpreter','Latex','Fontsize',20)
title('Iteration '+plotnam+', $\Delta=$ '+string(round(dist_curr_for+dist_curr_back,3)),'Interpreter','Latex','Fontsize',20)

subplot(1,2,2)
plot(PS_traj(:,num_spec + 1),PS_traj(:,num_spec + 2),'c-.','LineWidth',2)
hold on
plot(Ham_traj_for(:,num_spec+1),Ham_traj_for(:,num_spec+2),'LineWidth',2)
plot(Ham_traj_back(:,num_spec + 1),Ham_traj_back(:,num_spec+2),'b','LineWidth',2)
plot(PS_traj(:,num_spec + 1),PS_traj(:,num_spec + 3:end),'c-.','LineWidth',2)
plot(Ham_traj_for(:,num_spec+1),Ham_traj_for(:,num_spec+3:end),'LineWidth',2)
plot(Ham_traj_back(:,num_spec + 1),Ham_traj_back(:,num_spec+3:end),'b','LineWidth',2)
hold off
%grid on
h1 = legend('Descender','Fwd HamEoM','Bwd HamEoM','location','best');
set(h1,'Interpreter','latex');
h1.FontSize = 15;

%title('Descender vs HamEoM')
xlabel('$p_1$','Interpreter','Latex','Fontsize',20)
ylabel('$p_2$ onwards','Interpreter','Latex','Fontsize',20)

x0=0.1;
y0=0.1;
width=28;
height=8;
set(gcf,'position',[x0,y0,width,height])
set(gcf,'units','centimeters','position',[x0,y0,width,height])
if save_plot_flag == 1
    %save_plot_name = save_plot_name + '_' + string(pt_idx)+'.png';
    saveas(gcf,save_plot_name)
end

end