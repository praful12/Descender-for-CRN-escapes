function [] = plot_PS_traj_fft(PS_traj,t_traj,plotnam,save_plot_name,save_plot_flag)

    num_spec = size(PS_traj,2)/2;
    figure()
    subplot(1,2,1)
    hold on
    for i = 1:num_spec
        plot(abs(fft(PS_traj(:,i))),'LineWidth',2)
    end
    xlim([0 20])
    hold off
    grid on
    title(plotnam +': |fft(position)|')

    subplot(1,2,2)
    hold on
    for i = 1:num_spec
        plot(abs(fft(PS_traj(:,i+num_spec))),'LineWidth',2)
    end
    hold off
    xlim([0 20])
    grid on
    title(' |fft(momentum)|')

    x0=0.1;
    y0=0.1;
    width=23;
    height=8;

    set(gcf,'position',[x0,y0,width,height])
    set(gcf,'units','centimeters','position',[x0,y0,width,height])

    if save_plot_flag == 1
        saveas(gcf,save_plot_name)
    end

end