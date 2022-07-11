function [] = plot_PS_traj(PS_traj,t_traj,plotnam,save_plot_name,save_plot_flag)

    num_spec = size(PS_traj,2)/2;
    figure()
    subplot(1,3,1)
    plot(PS_traj(:,1),PS_traj(:,2:num_spec),'LineWidth',2)
    hold on
    plot(PS_traj(1,1),PS_traj(1,2:num_spec),'k*','MarkerSize',15)
    plot(PS_traj(end,1),PS_traj(end,2:num_spec),'bo','MarkerSize',15)
    hold off
    grid on
    xlabel('q_1')
    ylabel('q_2')
    title(plotnam +': Position')

    subplot(1,3,2)
    plot(PS_traj(:,1+num_spec),PS_traj(:,num_spec+2:2*num_spec),'LineWidth',2)
    hold on
    plot(PS_traj(1,1+num_spec),PS_traj(1,num_spec+2:2*num_spec),'k*','MarkerSize',15)
    plot(PS_traj(end,1+num_spec),PS_traj(end,num_spec+2:2*num_spec),'bo','MarkerSize',15)
    hold off
    grid on
    xlabel('p_1')
    ylabel('p_2')
    title('Momentum')
    
    subplot(1,3,3)
    plot(PS_traj(:,1),t_traj,'.-')
    hold on
    plot(PS_traj(1,1),t_traj(1),'k*','MarkerSize',15)
    plot(PS_traj(end,1),t_traj(end),'bo','MarkerSize',15)
    hold off
    grid on
    xlabel('q_1')
    ylabel('Time')
    title('Time along trajectory')
    
    x0=0.1;
    y0=0.1;
    width=28;
    height=8;
    set(gcf,'position',[x0,y0,width,height])
    set(gcf,'units','centimeters','position',[x0,y0,width,height])

    if save_plot_flag == 1
        saveas(gcf,save_plot_name)
    end

end