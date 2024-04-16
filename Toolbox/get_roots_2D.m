%Get roots using mass action kinetics for 2D models. This will not work for any high
%dimensianal problem. For that use Bertini and get roots.
%Plot Roots on morse-function
%Input: MAK
%Output: eig_val_lst, pos_root_arr, unst_vec, stab_idx, sad_idx, x_ic_arr

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
%pos_root_arr

%Get eigenvalues of jacobian and assign stability
eig_val_lst = [];
c = sym('c',[1,num_spec]);
stab_idx = [];
sad_idx = [];
unst_vec = [];
x_ic_arr = [];  %For getting initial conditions
eps = 10^(-3);

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

    if pos_eval_ct == 0    %stable root
        stab_idx = [stab_idx ; i];
    elseif pos_eval_ct == 1 %saddle root
        sad_idx = [sad_idx ; i];
        sad_root = pos_root_arr(i,:);
        x_ic_arr = [ x_ic_arr; sad_root + eps*unst_vec(end)];
        x_ic_arr = [ x_ic_arr; sad_root - eps*unst_vec(end)];
    end
end


eig_val_lst;
new_x_ic_arr = zeros(size(x_ic_arr));
new_x_ic_arr(:,:)=x_ic_arr;
x_ic_arr = new_x_ic_arr;
pos_root_arr
x_ic_arr

%{
figure()
%lb = ; ub = 4;
hold on
plot(pos_root_arr(stab_idx(1),1),pos_root_arr(stab_idx(1),2),'k*','MarkerSize',20,'MarkerFaceColor','k')
plot(pos_root_arr(sad_idx(1),1),pos_root_arr(sad_idx(1),2),'bo','MarkerSize',10,'MarkerFaceColor','b')
plot(pos_root_arr(2,1),pos_root_arr(2,2),'rs','MarkerSize',10,'MarkerFaceColor','r')
for i = 2:size(stab_idx,1)
    plot(pos_root_arr(stab_idx(i),1),pos_root_arr(stab_idx(i),2),'k*','MarkerSize',20)
end
for i = 2:size(sad_idx,1)
    plot(pos_root_arr(sad_idx(i),1),pos_root_arr(sad_idx(i),2),'bo','MarkerSize',10,'MarkerFaceColor','b')
end
hold off
%legend('Stable roots','Saddle roots','All roots','location','north')
h1=legend('Stable','Saddle','Bi-saddle','location','best');
%h1=legend('Stable','Saddle','location','best');
set(h1,'Interpreter','latex');
h1.FontSize = 15;
ax=gca;
ax.FontSize = 13;

xlabel('$q_1$','interpreter','latex','FontSize',20)
ylabel('$q_2$','interpreter','latex','FontSize',20)
lim = double([min(pos_root_arr(:,1))-0.2  max(pos_root_arr(:,1))+0.2]);
xlim(lim)
lim = double([min(pos_root_arr(:,2))-0.2  max(pos_root_arr(:,2))+1]);
%ylim(lim)
%grid on
title('2-Schlogl : Stationary points','Interpreter','latex','FontSize',20)

%}

%Morse function
%{

x = linspace(lb,ub);
y = linspace(lb,ub);
mf = int(MAK(1),q(1))+int(MAK(2),q(2));
mf_func =  matlabFunction(mf,'vars',{q});
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));

for i = 1:size(Z,1)
for j = 1:size(Z,2)
Z(i,j) = mf_func([X(i,j) Y(i,j)]);
end
end

contour(X,Y,Z,40)
hold on
for i = 1:size(stab_idx,1)
    plot(pos_root_arr(stab_idx(i),1),pos_root_arr(stab_idx(i),2),'k*','MarkerSize',20)
end
for i = 1:size(sad_idx,1)
    plot(pos_root_arr(sad_idx(i),1),pos_root_arr(sad_idx(i),2),'bo','MarkerSize',20)
end
hold off


for i = 1:size(Z,1)
for j = 1:size(Z,2)
Z(i,j) = subs(mf,q,[x(i) y(j)]);
end
end
contour(X,Y,Z)

%}
%{
figure()
%fsurf(int(poly_ode(1),c(1))+int(poly_ode(2),c(2)),[40 90 0 30])

contour(-int(MAK(1),q(1))-int(MAK(2),q(2)))
%}
%{
hold on
for i = 1:size(stab_idx,1)
    plot(pos_root_arr(stab_idx(i),1),pos_root_arr(stab_idx(i),2),'*','MarkerSize',20)
end
for i = 1:size(sad_idx,1)
    plot(pos_root_arr(sad_idx(i),1),pos_root_arr(sad_idx(i),2),'o')
end
hold off
shading(gca,'interp')
%set(gca,'xtick',[1:num_spec],'xticklabel',species_uq)
xlim([lb ub])
ylim([lb ub])
%ylabel('Concentration')
%xlabel('Species')
grid on
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
