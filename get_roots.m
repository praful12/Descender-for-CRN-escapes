%Get roots using mass action kinetics. This will not work for any high
%dimensianal problem. For that use Bertini and get roots.

roots = vpasolve(MAK,q,ones(1,num_spec));
root_arr = subs(q,roots);
num_root = size(root_arr,1);
%
%find positive roots
pos_root_arr = [];
ct = 0;
for i = 1:num_root
    for j = 1:num_spec
        if root_arr(i,j)>0
            ct = ct + 1;
        end
    end
    if ct == num_spec
        pos_root_arr = [pos_root_arr;root_arr(i,:)];
    end
    ct = 0;
end
num_pos_root = size(pos_root_arr,1);
pos_root_arr

%Get eigenvalues of jacobian and assign stability
eig_val_lst = [];
c = sym('c',[1,num_spec]);
stab_idx = [];
sad_idx = [];
unst_vec = [];

for i = 1:num_pos_root
    pos_eval_ct = 0;
    [V,D] = eig(jac_fun(pos_root_arr(i,:)));
    eig_val_lst = [eig_val_lst; double(diag(D)')];
    for j = 1:num_spec
        if real(D(j,j))>0
            pos_eval_ct = pos_eval_ct + 1;
            unst_vec = [unst_vec ; V(:,j)'];
        end
    end
    if pos_eval_ct == 0
        stab_idx = [stab_idx ; i];
    elseif pos_eval_ct == 1
        sad_idx = [sad_idx ; i];
    end
end
eig_val_lst;
pos_root_arr;
%}

%
%plot stable and saddle roots
figure()
plot([1:num_spec],pos_root_arr(stab_idx(1),:),'k*-','MarkerSize',10)
hold on
plot([1:num_spec],pos_root_arr(sad_idx(1),:),'bo-','MarkerSize',10)
%plot(pos_root_arr(:,1),pos_root_arr(:,2),'r.','MarkerSize',10)
for i = 2:size(stab_idx,1)
    plot([1:num_spec],pos_root_arr(stab_idx(i),:),'k*-','MarkerSize',10)
end
for i = 2:size(sad_idx,1)
    plot([1:num_spec],pos_root_arr(sad_idx(i),:),'bo-','MarkerSize',10)
end
hold off
legend('Stable roots','Saddle roots','location','best')
set(gca,'xtick',[1:num_spec],'xticklabel',species_uq)
xlim([0 num_spec+1])
ylabel('Concentration')
xlabel('Species')
grid on
title(plotnam)
%title("Stable(circle) and saddle(square) "+string(model_name),'Interpreter','none')
saveas(gcf,'../4 Plot and publish/'+model_name+'_stability.png')
%}

%{
%plot positive roots
figure()
hold on
for i = 1:num_pos_root
    plot([1:num_spec],pos_root_arr(i,:),'o-')
end
hold off
set(gca,'xtick',[1:num_spec],'xticklabel',species_uq)
xlim([0 num_spec+1])
ylabel('Concentration')
xlabel('Species')
grid on
title("Positive roots of "+string(model_name) + " : total #"+string(num_pos_root),'Interpreter','none')
saveas(gcf,'../4 Plot and publish/'+model_name+'_pos_roots.png')
%}
