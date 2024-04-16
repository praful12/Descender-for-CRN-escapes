function flag = verify_const(iter, P_val,xtraj,Ham_fun)
flag = 0;
if P_val(1) < 0
    "TIme negative - constraints not satisfied in idx "+string(iter)
    flag = 1;
elseif abs(Ham_fun([xtraj,P_val(2:end)]))>10^(-3)
    "Hamiltonian ~= 0 in idx " + string(iter)
    flag = 1;
end
    
end