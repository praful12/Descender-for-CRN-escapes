
%Selkov model
num_spec = 2;
comp_vec = [[0,0];[1,0];[1,2];[0,3];[0,1]];

alpha = -1/3;
k2 = (1-alpha^2)/(2*14^2); k2_b = k2;
k3 = 1- alpha; k3_b = 2*(1-alpha)*(36/49);
k1 = 4*14-k3_b; k1_b = 1 + alpha;

rxn_rate = [[1,2,k1];[2,1,k1_b];[3,4,k2];[4,3,k2_b];...
    [5,1,k3];[1,5,k3_b]];

%rxn_rate
%------------------------------------------------------------
%Hypergraph
comp_vec_name = string(zeros(size(comp_vec,1),1));
comp_vec_name = string(string(comp_vec(:,1))+ 'X_'+string(1));
for i = 2:num_spec
    comp_vec_name = comp_vec_name + '+' +   string(comp_vec(:,i))+ 'X_'+string(i) ;
end
comp_vec_name;


species_uq = [];
for i = 1:num_spec
    species_uq = [species_uq , 'X_'+string(i)];
end
species_uq;