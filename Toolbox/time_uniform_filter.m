%Filter the trajectory in time uniform sampling

function traj_n_s = time_uniform_filter(traj,t_traj,traj_pts_new,f_filter)

    t_uniform = t_traj(1):(t_traj(end)-t_traj(1))/(traj_pts_new-1):t_traj(end); t_uniform = t_uniform';
    %Time uniform sample
    traj_n = interp1(t_traj,traj,t_uniform);
    s_ran = 0:1/(traj_pts_new-1):1;
    traj_st = traj_n(end,:) + (1-s_ran').*(traj_n(1,:)-traj_n(end,:));
    traj_left = traj_n - traj_st;
    
    %Smoothen difference from straight line and add back to get a smoothened
    %trajectory
    traj_left_s = smooth_traj_butter(traj_left,f_filter);
    traj_n_s = traj_st + traj_left_s;
    
    %{
    figure()
    subplot(1,1,1)
    plot(traj_n_s(:,1),traj_n_s(:,2:end))
    hold on
    plot(traj(:,1),traj(:,2:end),'--')
    hold off
    legend('filtered','original','location','best')
    
    subplot(1,2,2)
    semilogy(abs(fft(traj_n_s(:,1))).^2)
    hold on
    semilogy(abs(fft(traj(:,1))).^2,'--')
    semilogy(abs(fft(traj_n_s(:,2:end))).^2)
    semilogy(abs(fft(traj(:,2:end))).^2,'--')
    hold off
    legend('filtered','original','location','best')
    xlim([0 20])
    title('Spectrum')

    x0=0.1;
    y0=0.1;
    width=28;
    height=8;
    set(gcf,'position',[x0,y0,width,height])
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
    %}
    
end