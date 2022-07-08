function P_val = eval_mom_num(xtraj,dxtraj,eqns_fsolve_fun,Jac_eqns_fsolve_fun)
    
    num_spec = size(xtraj,2);
    %---------- FSOLVE/LSQNONLIN inputs ------    
    %initial guess for minimizer
    var_ic = 0.1*ones(1,num_spec+1);
    %Bounds
    lb = -100* ones(1,num_spec+1);
    lb(1) = 0;
    ub = 100*ones(1,num_spec+1);
    ub(1) = 1000;


    %solve function and gradient 
    fun_fsolve_handle = @(X) fun_fsolve([X(1),xtraj,X(2:num_spec+1),dxtraj],eqns_fsolve_fun,Jac_eqns_fsolve_fun);
    
    %Solve function
    options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off','TolFun',10^(-12),'MaxIter',10^8,'MaxFunEvals',10^8);
    %[P_val,fval,exitflag,output] = fsolve(fun_fsolve_handle,var_ic,options);
    [P_val,~,~,~] = lsqnonlin(fun_fsolve_handle,var_ic,lb,ub,options);


end