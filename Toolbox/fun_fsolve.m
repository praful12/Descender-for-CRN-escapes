function [f_vec, f_jac] = fun_fsolve(X,eqns,Jac_eqns)
    f_vec = eqns(X);
    f_jac = Jac_eqns(X);
end