function [MAK, jac, MAK_fun, jac_fun] = get_MAK(dHamdp,Ham_qp_mix,num_spec)

p = sym('p', [1 num_spec]);
q = sym('q' ,[1 num_spec]);

%--- get polynomails for Mass-Action kinetics and its Jacobian
MAK = subs(dHamdp,p,zeros(1,num_spec));
jac = subs(Ham_qp_mix,p,zeros(1,num_spec));

%---------Extract MATLAB function handles-------
MAK_fun = matlabFunction(MAK,'vars',{q});
jac_fun = matlabFunction(jac,'vars',{q});
