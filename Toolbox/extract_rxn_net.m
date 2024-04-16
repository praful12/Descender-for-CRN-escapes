function [num_spec,species,comp_vec,rxn_rate,G] = extract_rxn_net(filnam)

%File name
fileID = fopen(filnam + ".txt");

%read file
rxn_net = {};
line = fgetl(fileID);

%count empty lines - this is if there are many reaction networks in a file
%separated by an empty line
emp_ct = 0;

while ischar(line)
    line = strtrim(line);
    if size(line,2) == 0 
        emp_ct = emp_ct + 1;
        if emp_ct == 1
            break
        end
    end
    rxn_net{end+1} = line;
    
    line = fgetl(fileID);
end
rxn_net;


%get input and output complex from each line
in_comp = {};
out_comp = {};
rate_c_val = {};

for i = 1 : size(rxn_net,2)
    %input complex
    sp_str = split(rxn_net{i},'-');
    in_comp{end+1} = strtrim(sp_str{1});
    
    %rate constant
    ct = 0;
    for j = 1:size(sp_str,1)
        if size(sp_str{j}) ~= 0
            ct = ct + 1;
            if ct > 1 && size(sp_str{j},1)>0
                rate_c_val{end+1} = sp_str{j};
                if size(sp_str{j+1},1)>0
                    mp = strsplit(string(sp_str(j+1)));
                    if mp(1)~=">" && size(mp,1)>0
                        rate_c_val{end} = string(rate_c_val{end}+"-"+string(sp_str(j+1)));
                    end
                end
                break
            end
        end
    end
    rate_c_val{end} = str2double(erase(rate_c_val{end}," "));
    %out complex
    sp_str = split(rxn_net{i},'>');
    out_comp{end+1} = strtrim(sp_str{end});
    
    
end

%get a graph of rxn_network
G = digraph(in_comp,out_comp);

%Extract species and stoichiometry from nodes of the graph or complexes
%ascii 0 = 48 A = 65
G_nodes = G.Nodes;
stoich = {};
species = {};
comp_sz = [];

for i = 1:size(G_nodes,1)
     str_nam = strtrim(split(G_nodes{i,1},'+'));
     str_sp = char(str_nam);
     comp_sz = [comp_sz;size(str_nam,1)];
     for j = 1:size(str_nam,1)
         dl = 0;
         for k = 1:size(str_sp(j,:),2)
             if str_sp == '0'
                 stoich{end+1} = '0';
                 species{end+1} = '0';
                 break
             end
             if double(str_sp(j,k))>64
                dl = k;
                species{end+1} = str_sp(j,k:end);
                if dl == 1
                    stoich{end+1} = '1';
                else
                     stoich{end+1} = str_sp(j,1:dl-1);
                end
                break
             end
         end
         species{end} = strtrim(species{end});
     end
end
stoich;
species;
%

%get unique species. c1 is names, c2 is index
[species_uq, c1 ,c2] = unique(species);
num_spec = size(species_uq,2);
cum_comp = cumsum(comp_sz);
cum_comp = [0;cum_comp];


%now write complexes as vectors
comp_vec = [];
flag = 1;
for i = 1:size(comp_sz,1)
    vec = zeros(1,num_spec);
    for j = 1:comp_sz(i)
        sp_1 = sscanf(stoich{1,cum_comp(i)+j},'%i');
%        sp_2 = sscanf(species{1,j+cum_comp(i)},'%s');
        vec(1,c2(j+cum_comp(i))) = sp_1;
    end    
    comp_vec = [comp_vec;vec];
end
comp_vec;

%Get reaction matrix to go into Extract_Hamiltonian
num_rxn = size(rate_c_val,2);
rxn_rate = zeros(num_rxn,3);
for i = 1:num_rxn
    rxn_rate(i,1) = find(strcmp(in_comp{i},G_nodes{:,1}));
    rxn_rate(i,2) = find(strcmp(out_comp{i},G_nodes{:,1}));
    rxn_rate(i,3) = rate_c_val{1,i};
end

rxn_rate;

%remove zero complex
if species_uq{1} =='0'
    num_spec = num_spec - 1;
    comp_vec = comp_vec(:,2:end);
    species_uq(1) = [];
end