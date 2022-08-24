%Function to plot unsmoothed and smoothed delta_x
%set save_plot_flag=1 for saving the plot and give save_plot_name
function c_out = debugging_deltax_plots(delta_x,delta_x_s,f0,plotnam,save_plot_name,save_plot_flag)
    num_spec = size(delta_x,2);
    figure()
    subplot(1,2,1)
    hold on
    plot(delta_x_s(:,1),'-','Linewidth',2)
    plot(delta_x(:,1),'--','Linewidth',1.8)
    for i = 2:num_spec
        plot(delta_x_s(:,i),'-','Linewidth',2)
        plot(delta_x(:,i),'--','Linewidth',1.8)
    end
    %legend('smoothened','original','location','best')
    ax=gca;
    ax.FontSize = 15;
    xlabel('Sample point along curve','interpreter','latex','FontSize',20)
    ylabel('$g$ and $g_s$','interpreter','latex','FontSize',20)
    title('Functional gradient','interpreter','latex','FontSize',20)
    hold off

    subplot(1,2,2)
    semilogy(abs(fft(delta_x_s(:,1))).^2,'-','Linewidth',2)
    hold on
    plot(abs(fft(delta_x(:,1))).^2,'--','Linewidth',1.8)
    for i = 2:num_spec
        plot(abs(fft(delta_x_s(:,i))).^2,'-','Linewidth',2)
        plot(abs(fft(delta_x(:,i))).^2,'--','Linewidth',1.8)
    end
    %grid on
    h1 = legend('Filtered','Original','location','best');
    set(h1,'Interpreter','latex');
    h1.FontSize = 15;
    
    xlabel('Frequency','interpreter','latex','FontSize',20)
    ylabel('Power Spectrum','interpreter','latex','FontSize',20)
    title('Cutoff frequency : $f_c =$ '+string(f0),'interpreter','latex','FontSize',20)

    xlim([0 8*f0])
    %ylim([-0.5 0.5])
    hold off
    
    x0=0.1;
    y0=0.1;
    width=28;
    height=8;
    set(gcf,'position',[x0,y0,width,height])
    set(gcf,'units','centimeters','position',[x0,y0,width,height])

    if save_plot_flag == 1
        saveas(gcf,save_plot_name)
    end

    c_out = [];
end