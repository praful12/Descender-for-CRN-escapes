function PS_traj_fin = AFGD_4_hcnet(model_name,sad_pt,stab_pt,traj_pts,iter_max)

%load Hamiltonian
load('..\Data\'+model_name+'_Ham.mat')

%load initial condition
load('..\Data\' + model_name + '_initcond_'+string(sad_pt)...
        +'_'+string(stab_pt)+'straight.mat')
traj_ic = double(traj_ic);
[traj,len_curve_prime] = space_uniform_sampling(traj_ic,traj_pts);


%Initialize arrays for saving and delta_x_s
Initialize_descender

%Save to file
savenam = '..\Data\' + model_name + '_hc_network_'+string(sad_pt)...
+'_'+string(stab_pt)+'straight_'+string(traj_pts)+'_'+string(date)+'.mat';


%Actual descender
eps_ic = 1; eps = eps_ic; eps_min = 0.1;     %Initial and least step size
delta_S = 10; err_thresh = 10^(-1); 
f0_max = 40; f0_step = 0.1; f0_init = f0;
eps_S_thresh = 10^(-7); delta_S_thresh = 5*10^(-7);
iter = 1; a_min=100; a_thresh = 0.001;
delta_S_thresh_eps = 10^(-7); time_uniform_f = 45;


%Descend
while f0<f0_max  &&  iter<iter_max % && a_min>a_thresh
    iter = iter + 1;
    traj_pts = size(traj,1);

    %take a step in the gradient by picking appropriate step size
    [eps, action_S_new, PS_traj, t_traj, delta_x, S_traj] = eps_picker...
        (traj,eqns_fsolve_fun,Jac_eqns_fsolve_fun,dHamdq_fun,Ham_fun_qp,delta_x_s,action_S,eps,num_smooth,eps_min,delta_S_thresh_eps);
    traj = PS_traj(:,1:num_spec);
    delta_S = action_S_new - action_S;
    [iter, eps, action_S_new, delta_S, f0]
    eps_used = eps;

    

    if abs(delta_S)<delta_S_thresh || eps < eps_min
        f0 = f0 + 0.1; eps = eps_ic; delta_S = 1;  
        traj = time_uniform_filter(traj,t_traj,traj_pts,time_uniform_f);
        [traj,len_curve_prime] = space_uniform_sampling(traj,traj_pts);
        [PS_traj, t_traj, delta_x] = traj_2_PS_traj_func(traj,eqns_fsolve_fun,Jac_eqns_fsolve_fun,dHamdq_fun,Ham_fun_qp,num_smooth);    
    end

    
    %smoothen gradient for next iteration
    delta_x_s =  smooth_traj_butter(delta_x,f0); 

    %Find closest distance by integrating HamEoM
    %a_min = Ham_closest_approach_noplot(PS_traj,S_traj,t_traj,dHamdp_fun, dHamdq_fun,save_plot_name,1,plotnam)

    
    %Update descender and save to file
    Update_descender
    save(savenam)
    

    
end
PS_traj_fin = PS_traj;
end
