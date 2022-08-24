%Plot the heteroclinic network trajecotries simultaneously on 2 indices
%Also plot escape and delta_G s

%{
idx1 = 1; idx2 = 2;
x_y_idx = size(hc_traj_arr,3)*0.4;
figure()
hold on
for i = 1:num_sad
    for j = 1:2
        plot(squeeze(hc_traj_arr(i,j,:,idx1)),squeeze(hc_traj_arr(i,j,:,idx2)),'--')
        plot(squeeze(hc_esc_traj_arr(i,j,:,idx1)),squeeze(hc_esc_traj_arr(i,j,:,idx2)),'-')
        in_p = squeeze(hc_traj_arr(i,j,1,:));
        out_p = squeeze(hc_traj_arr(i,j,end,:));
        mid_p = (in_p + out_p)/2;
        text(mid_p(1),mid_p(2),"\Deltat_e="+string(round(hc_esc_time(i,j))))

        %mid_p = min(in_p,out_p);
        %text(mid_p(1),mid_p(2),"\Deltat="+string(hc_esc_time(i,j))+";\DeltaG="+string(hc_delta_G(i,j)))
        plot(squeeze(hc_traj_arr(i,j,1,idx1)),squeeze(hc_traj_arr(i,j,1,idx2)),'bo','MarkerSize',20)
        plot(squeeze(hc_traj_arr(i,j,end,idx1)),squeeze(hc_traj_arr(i,j,end,idx2)),'k*','MarkerSize',20)
    end
    xlim([1 3])
    ylim([1 3])
    legend('Relaxation','Escape','Saddle','Stable','location','best')
    %grid on
    title("2-Schlogl : kdif = "+string(kdif))
end
hold off
saveas(gcf,'../4 Plot and publish/'+model_name+'_esc_hcnet.png') 
%}
%{
N=4;
figure()
hold on
for i = 1:num_sad
    for j = 1:2
        traj =     squeeze(hc_esc_traj_arr(i,j,:,:));
        plot(squeeze(hc_esc_traj_arr(i,j,:,idx1)),squeeze(hc_esc_traj_arr(i,j,:,idx2)),'-')
        c1 = squeeze(hc_traj_arr(i,j,end,idx1))+0.5*(squeeze(hc_traj_arr(i,j,1,idx1))-squeeze(hc_traj_arr(i,j,end,idx1)));
        c2 = squeeze(hc_traj_arr(i,j,end,idx2))+0.5*(squeeze(hc_traj_arr(i,j,1,idx2))-squeeze(hc_traj_arr(i,j,end,idx2)));
        text(c1,c2,"\DeltaG="+string( round(hc_delta_G(i,j),N)))
        plot(squeeze(hc_traj_arr(i,j,1,idx1)),squeeze(hc_traj_arr(i,j,1,idx2)),'bo','MarkerSize',12)
        plot(squeeze(hc_traj_arr(i,j,end,idx1)),squeeze(hc_traj_arr(i,j,end,idx2)),'k*','MarkerSize',15)
    end
    ylim([1 3.2])
    xlim([1 3.2])
    legend('Escape','Saddle','Stable','location','best')
    %grid on
    title("2-Schlogl : kdif = "+string(kdif))
end
hold off
saveas(gcf,'../4 Plot and publish/'+model_name+'_deltaG.png')
%}
%{
N=4;
S_traj_arr = zeros(num_sad,2,traj_pts);
for i = 1:num_sad
    for j = 1:2
        traj =     squeeze(hc_esc_traj_arr(i,j,:,:));
        [PS_traj, t_traj, delta_x] = traj_2_PS_traj_func(traj,err_thresh,eqns_fsolve_fun,Jac_eqns_fsolve_fun,dHamdq_fun);
        [action_S_new, S_traj] = action_PS_traj(PS_traj);
        S_traj_arr(i,j,:) = S_traj;
    end
end


S_traj_arr_or = S_traj_arr;
curr_sad = 1;
curr_stab = 2;
val_at_stab = 0; 
hc_net;

for i = 1:num_sad
    S_traj_arr(curr_sad,curr_stab,:) =  S_traj_arr(curr_sad,curr_stab,:) + val_at_stab - S_traj_arr(curr_sad,curr_stab,1);

 %   [hc_net(curr_sad,1) hc_net(curr_sad,curr_stab+1)]
 %   [S_traj_arr(curr_sad,curr_stab,1) S_traj_arr(curr_sad,curr_stab,end)]
    
    val_at_sad = S_traj_arr(curr_sad,curr_stab,end);
    curr_stab = mod(curr_stab,2) + 1;
    S_traj_arr(curr_sad,curr_stab,:) =  S_traj_arr(curr_sad,curr_stab,:) + val_at_sad - S_traj_arr(curr_sad,curr_stab,end);    
    
%    [hc_net(curr_sad,1) hc_net(curr_sad,curr_stab+1)]
%    [S_traj_arr(curr_sad,curr_stab,1) S_traj_arr(curr_sad,curr_stab,end)]
    val_at_stab = S_traj_arr(curr_sad,curr_stab,1);
  
    for j = 1:num_sad
        if j ~= curr_sad && min(abs(hc_net(curr_sad,curr_stab+1)*[1 1 1]-hc_net(j,:)))==0
                %hc_net(j,:)
                %abs(hc_net(curr_sad,curr_stab+1)*[1 1 1]-hc_net(j,:))
                [~,curr_stab] = min(abs(hc_net(curr_sad,curr_stab+1)*[1 1 1]-hc_net(j,:)));
                curr_stab = curr_stab-1;
                curr_sad = j ;              
                %[curr_sad curr_stab]
                break
                
        end
    end

end
%}

figure()
for i = 1:num_sad
    for j = 1:2
        traj = squeeze(hc_esc_traj_arr(i,j,:,:));
        S_traj = squeeze(S_traj_arr(i,j,:));
        plot3(traj(:,1),traj(:,2),S_traj,'Linewidth',2)
        hold on
        plot(squeeze(hc_esc_traj_arr(i,j,:,idx1)),squeeze(hc_esc_traj_arr(i,j,:,idx2)),'--','Linewidth',2)
        c1 = squeeze(hc_traj_arr(i,j,end,idx1))+0.5*(squeeze(hc_traj_arr(i,j,1,idx1))-squeeze(hc_traj_arr(i,j,end,idx1)));
        c2 = squeeze(hc_traj_arr(i,j,end,idx2))+0.5*(squeeze(hc_traj_arr(i,j,1,idx2))-squeeze(hc_traj_arr(i,j,end,idx2)));
%        text(c1,c2,"\DeltaG="+string( round(hc_delta_G(i,j),N)))
        plot(squeeze(hc_traj_arr(i,j,1,idx1)),squeeze(hc_traj_arr(i,j,1,idx2)),'bo','MarkerSize',12)
        plot(squeeze(hc_traj_arr(i,j,end,idx1)),squeeze(hc_traj_arr(i,j,end,idx2)),'k*','MarkerSize',15)
    end
    ylim([1 3.2])
    xlim([1 3.2])
    legend('Action','Escape','Saddle','Stable','location','best')
    view([-66.9000 19.2000])
    grid on
    xlabel('q_'+string(idx1))
    ylabel('q_'+string(idx2))
    zlabel('Action')
    title("2-Schlogl : kdif = "+string(kdif))
end
hold off

saveas(gcf,'../4 Plot and publish/'+model_name+'_deltaG_3D.png')
%}