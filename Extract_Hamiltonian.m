
%-----------------------Hamiltonian ----------------------------
p = sym('p', [1 num_spec]);
q = sym('q' ,[1 num_spec]);

Ham = sym(0);

for i = 1:size(rxn_rate,1)
    y_in = comp_vec(rxn_rate(i,1),:)';
    y_out = comp_vec(rxn_rate(i,2),:)';
    Ham = Ham + rxn_rate(i,3)*(exp(p*(y_out-y_in))-1)*prod(q.'.^y_in);
end
Ham

%------- Gradient ------- 
dHamdp = gradient(Ham,p);
dHamdq = gradient(Ham,q);

%--- get ode for Mass-Action kinetics
MAK = subs(dHamdp,p,zeros(1,num_spec));

%-------------------- higher order derivatives ----------------
Ham_q_hessian = jacobian(dHamdq,q);
Ham_p_hessian = jacobian(dHamdp,p);
Ham_qp_mix = jacobian(dHamdq,p);

jac = subs(Ham_qp_mix,p,zeros(1,num_spec));

%------- Non-linear optimization ------- 
%Hamiltonian is convex in momentum. 

syms dx [1 num_spec]
syms dt
coord_min = [dt p];

%Minimization function, gradient and hessian
func_min = Ham*dt - p*dx.' ;
func_min_grad = gradient(func_min,coord_min);
func_min_hess = jacobian(func_min_grad,coord_min);

%Constraint function, gradient and hessian
func_con = Ham;
func_con_grad = gradient(func_con,coord_min);
func_con_hess = jacobian(func_con_grad,coord_min);

%-------------- fsolve equations --------------  
eqns_fsolve = [Ham; dHamdp*dt - dx.' ];
Jac_eqns_fsolve = jacobian(eqns_fsolve,coord_min);

%---------Extract MATLAB function handles-------
%For Heteroclinic network
MAK_fun = matlabFunction(MAK,'vars',{q});
jac_fun = matlabFunction(jac,'vars',{q});

%For descent
PS_coord = [q p];
Ham_fun_qp = matlabFunction(Ham,'vars',{PS_coord});
dHamdq_fun = matlabFunction(dHamdq','vars',{PS_coord});
dHamdp_fun = matlabFunction(dHamdp','vars',{PS_coord});
Ham_q_hess_fun = matlabFunction(Ham_q_hessian,'vars',{PS_coord});
Ham_p_hess_fun = matlabFunction(Ham_p_hessian,'vars',{PS_coord});
Ham_qp_mix_fun = matlabFunction(Ham_qp_mix,'vars',{PS_coord});

%For minimization
coord_prime = [dt q p dx];
min_fun =  matlabFunction(func_min,'vars',{coord_prime});
min_grad_fun = matlabFunction(func_min_grad,'vars',{coord_prime});
min_hess_fun = matlabFunction(func_min_hess,'vars',{coord_prime});
con_fun =  matlabFunction(func_con,'vars',{coord_prime});
con_grad_fun = matlabFunction(func_con_grad,'vars',{coord_prime});
con_hess_fun = matlabFunction(func_con_hess,'vars',{coord_prime});

%For fsolve
eqns_fsolve_fun = matlabFunction(eqns_fsolve,'vars',{coord_prime});
Jac_eqns_fsolve_fun = matlabFunction(Jac_eqns_fsolve,'vars',{coord_prime});


%---------- Get equations for analytic evaluation -----
p = sym('p', [1 num_spec]);
q = sym('q' ,[1 num_spec]);
z=sym('z',[1,num_spec],['positive']);

Ham_qz = sym(0);

for i = 1:size(rxn_rate,1)
    y_in = comp_vec(rxn_rate(i,1),:);
    y_out = comp_vec(rxn_rate(i,2),:);
    Ham_qz = Ham_qz + rxn_rate(i,3)*(prod(z.^(y_out-y_in))-1)*prod(q.^y_in);
end

%------- Gradient ------- 
dHamdz = gradient(Ham_qz,z);
dq = sym('dq',[1,num_spec]);

Eqns = [Ham_qz==0];
for i = 2:num_spec
    Eqns = [Eqns; dHamdz(i)/dHamdz(1)*z(i)/z(1) == dq(i)/dq(1)];
end
Eqns


%{
%For conversion to coherent state
a=sym('a',[1,num_spec],['positive']);
q_NP = z.*a;
p_NP = log(z);
Ham_CS = subs(Ham,[p ; q],[p_NP;q_NP])
%}
