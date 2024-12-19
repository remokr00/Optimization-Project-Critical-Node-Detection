function out=ConnectedNP(G,selected)
% ConnectedNP(G,selected) - Computes the number of node pairs that can communicate
%         if the selected nodes are eliminated (returns -1 for invalid input data)
%
% G:         graph of the network
% selected:  a row array with IDs of selected nodes
    
    nNodes= numnodes(G);
    if length(selected)>=1
        if (max(selected)>nNodes || min(selected)<1 || length(unique(selected))<length(selected))
            out= -1;
            return
        end
    end
    aux= setdiff(1:nNodes,selected);
    Gr= subgraph(G,aux);
    dist= distances(Gr);
    out= (sum(dist(:)<Inf) - numnodes(Gr))/2;
end