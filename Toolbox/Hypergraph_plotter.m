%Hypergraph

G = digraph();
G = addnode(G,comp_vec_name);
for i = 1:size(rxn_rate,1)
    G = addedge(G,comp_vec_name(rxn_rate(i,1)),comp_vec_name(rxn_rate(i,2)),rxn_rate(i,3));
end
figure()
h = plot(G);
text(h.XData+.02, h.YData+.01 ,h.NodeLabel, ...
    'VerticalAlignment','Bottom',...
    'HorizontalAlignment', 'left',...
    'FontSize', 10)
h.NodeLabel = {};
title(plotnam+' reaction network')

%title(model_name + ' Hypergraph','Interpreter','none')