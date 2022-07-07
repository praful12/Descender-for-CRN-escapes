%Schlogl Model

%define complexes
num_spec = lattice_edge;
comp_vec = [zeros(1,num_spec)];
for i = 1:num_spec
    cv = zeros(1,num_spec);
    cv(i) = 1;
    comp_vec = [comp_vec; cv; 2*cv; 3*cv];
end
%comp_vec

%------------------------------------------------------------
%reaction network
rxn_rate = [];

%Schlogl for each site
for i = 1:num_spec
    c_in = 1;
    c_out = 1 + 3*(i-1) + 1;
    rate = r01;
    rxn_rate = [rxn_rate;[c_in c_out rate]];

    c_out = 1;
    c_in = 1 + 3*(i-1) + 1;
    rate = r10;
    rxn_rate = [rxn_rate;[c_in c_out rate]];
    
    c_in = 1 + 3*(i-1) + 2;
    c_out = 1 + 3*(i-1) + 3;
    rate = r23;
    rxn_rate = [rxn_rate;[c_in c_out rate]];
    
    c_in = 1 + 3*(i-1) + 3;
    c_out = 1 + 3*(i-1) + 2;
    rate = r32;
    rxn_rate = [rxn_rate;[c_in c_out rate]];    
end

%Diffusion
for i = 1:num_spec
    for j = i+1:num_spec
        c_in = 1 + 3*(i-1) + 1;
        c_out = 1 + 3*(j-1) + 1;
        rate = kdif;
        rxn_rate = [rxn_rate;[c_in c_out rate]];

        c_in = 1 + 3*(j-1) + 1;
        c_out = 1 + 3*(i-1) + 1;
        rate = kdif;
        rxn_rate = [rxn_rate;[c_in c_out rate]];

    end
end


%rxn_rate

%------------------------------------------------------------
%Names for hypergraph
comp_vec_name = string(string(comp_vec(:,1))+ 'X_'+string(1));
for i = 2:num_spec
    comp_vec_name = comp_vec_name + '+' +   string(comp_vec(:,i))+ 'X_'+string(i) ;
end
comp_vec_name;

%comp_vec_name

species_uq = [];
for i = 1:num_spec
    species_uq = [species_uq , 'X_'+string(i)];
end
species_uq;
