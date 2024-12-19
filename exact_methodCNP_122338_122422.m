%uploading the datas
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');
nNodes = size(Nodes,1); 
nLinks = size(Links,1);
G= graph(L);
c = [8, 10, 12];

%doing the operation for every c
for value = c
    %opening the file
    file_name = sprintf('exact_method_c_%d.lp', value);
    fid = fopen(file_name,'wt');
    %objective function
    fprintf(fid,'Min: ');
    for i=1:nNodes-1
        for j = i+1:nNodes
            if i ~= nNodes-1
                fprintf(fid,'+ u%d_%d ', i, j);
            else
                fprintf(fid,'+ u%d_%d;', i, j);
            end
        end
    end
    fprintf(fid,'\n');
    %1st costraint
    for i=1:nNodes
        fprintf(fid,'+ v%d ', i);
    end
    fprintf(fid,'= %d;\n',value);
    % 2nd constraint: Ensure that for each connected node pair, either one node is critical or the pair is considered disconnected
    for i = 1:size(L, 1)
        written = false;  % Initialize a flag to track if anything has been written for this 'i'
        for j = 1:size(L, 2)
            if L(i, j) ~= 0 && (i < j)  % Only consider each edge once (i < j)
                if ~written
                    written = true;  % Set flag to true as we are writing something
                    fprintf(fid, '+ u%d_%d + v%d + v%d >= 1;\n', i, j, i, j);
                else
                    fprintf(fid, '+ u%d_%d + v%d + v%d >= 1;\n', i, j, i, j);
                end
            end 
        end
    end
    % 3rd constraint: Ensure that for nodes not directly connected, their disconnection considers critical nodes and common neighbors
    for i = 1:nNodes
        for j = 1:nNodes
            if i ~= j && i < j && L(i, j) == 0  % Consider only non-connected pairs (i < j)
                ni = neighbors(G, i);  % Get neighbors of node i
                if ~isempty(ni)
                    for k = ni'  % Iterate over common neighbors
                        % Determine the correct form for u variables
                        if i < k
                            ui_k = sprintf('u%d_%d', i, k);
                        else
                            ui_k = sprintf('u%d_%d', k, i);
                        end
                        if j < k
                            uj_k = sprintf('u%d_%d', j, k);
                        else
                            uj_k = sprintf('u%d_%d', k, j);
                        end
                        if i < j
                            fprintf(fid, '+ u%d_%d - %s - %s - v%d >= - 1;\n', i, j, ui_k, uj_k, k);
                        else
                            fprintf(fid, '+ u%d_%d - %s - %s - v%d >= - 1;\n', j, i, ui_k, uj_k, k);
                        end
                    end
                end
            end
        end
    end
    
    %Binary variables costraints
    fprintf(fid, 'Bin ');
    for i=1:nNodes
        if i ~= nNodes
            fprintf(fid,'v%d,', i);
        else
            fprintf(fid,'v%d;', i);
        end
    end
    %closing the file
    fclose(fid);
end