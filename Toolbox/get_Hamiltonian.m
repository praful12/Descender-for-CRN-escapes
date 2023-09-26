function Ham = get_Hamiltonian(num_spec,comp_vec,rxn_rate)


%-----------------------Hamiltonian ----------------------------
p = sym('p', [1 num_spec]);
q = sym('q' ,[1 num_spec]);

Ham = sym(0);

for i = 1:size(rxn_rate,1)
    y_in = comp_vec(rxn_rate(i,1),:)';
    y_out = comp_vec(rxn_rate(i,2),:)';
    Ham = Ham + rxn_rate(i,3)*(exp(p*(y_out-y_in))-1)*prod(q.'.^y_in);
end

%Ham
