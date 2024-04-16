function [PS_traj, t_traj, delta_x] = plot_mom_deltax(traj,PS_traj,delta_x)

[traj_pts,num_spec] = size(traj);

figure()
plot(PS_traj(:,num_spec+1:end),'o-','Linewidth',1.2,'MarkerSize',5)
legend("p_1","p_2",'location','best')
title(plotnam + ': Momentum')
grid on
xlabel('Curve discretization index')
%saveas(gcf,'../4 Plot and publish/'+model_name+'_IC_momentum.png')

figure()
plot(delta_x(:,num_spec+1:end),'o-','Linewidth',1.2,'MarkerSize',5)
legend("delta-q_1","delta-q_2",'location','best')
title(plotnam + ': Momentum')
grid on
xlabel('Curve discretization index')


end
