function P_val = eval_mom_anal(xtraj,dxtraj,Eqns,dHamdp_fun)
    num_spec = size(xtraj,2);
    q = sym('q' ,[1 num_spec]);
    z=sym('z',[1,num_spec],['positive']);
    dq = sym('dq',[1,num_spec]);

    %-----------------Analytic VPASOLVE ----------------
    Eqn_2_solve = subs(Eqns,[q dq],[xtraj dxtraj]);
    P_val_solve = vpasolve(Eqn_2_solve,z);
    z_val = subs(z,P_val_solve);
    mom_val = log(z_val);

    %Pick the root along escape
    for j = 1:size(mom_val,1)
    vel = dHamdp_fun([xtraj mom_val(j,:)]);
    dt_val = dxtraj(1)/vel(1);
        if dt_val > 0 
            P_val = [dt_val mom_val(j,:)];
            break
        end
    end
    P_val = double(P_val);
end