% Importing the data
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');
nNodes = size(Nodes, 1);
nLinks = size(Links, 1);
G = graph(L);
time_limit = 60;
number_execution = 10;
c = [8, 10, 12]; % number of nodes in the solution

% GRASP parameter initialization
r = 5;
results_grasp = cell(length(c), 1);

% GA parameter initialization
population_size = 150;
q = 0.2; % mutation probability
num_elite = 15;

% Initialize arrays for stats and nodes corresponding to min and max values
min_grasp = zeros(size(c));
max_grasp = zeros(size(c));
mean_grasp = zeros(size(c));
min_nodes_grasp = cell(size(c));
max_nodes_grasp = cell(size(c));
iterations_grasp = cell(length(c), 1);

min_ga = zeros(size(c));
max_ga = zeros(size(c));
mean_ga = zeros(size(c));
min_nodes_ga = cell(size(c));
max_nodes_ga = cell(size(c));
iterations_ga = cell(length(c), 1);

% Main program to evaluate the GRASP algorithm
for idx = 1:length(c)
    results_grasp{idx} = zeros(number_execution, c(idx) + 2); % +2 to store the solution cost and iterations
    local_min_grasp = inf;
    local_max_grasp = -inf;
    iterations_grasp{idx} = zeros(number_execution, 1);

    % GRASP execution
    for j = 1:number_execution
        fprintf('Executing GRASP for c = %d. Iteration = %d\n', c(idx), j);
        [nodes, solution, iterations] = grasp(G, c(idx), r, time_limit);
        results_grasp{idx}(j, 1:c(idx)) = nodes;
        results_grasp{idx}(j, end-1) = solution;
        results_grasp{idx}(j, end) = iterations;
        iterations_grasp{idx}(j) = iterations;

        % Check for new min or max
        if solution < local_min_grasp
            local_min_grasp = solution;
            min_nodes_grasp{idx} = nodes;
        end
        if solution > local_max_grasp
            local_max_grasp = solution;
            max_nodes_grasp{idx} = nodes;
        end
    end

    % Compute min, mean, max
    min_grasp(idx) = local_min_grasp;
    mean_grasp(idx) = mean(results_grasp{idx}(:, end-1));
    max_grasp(idx) = local_max_grasp;
end

% Main program to evaluate the GA algorithm
for idx = 1:length(c)
    results_ga{idx} = zeros(number_execution, c(idx) + 2); % +2 to store the solution cost and iterations
    local_min_ga = inf;
    local_max_ga = -inf;
    iterations_ga{idx} = zeros(number_execution, 1);

    % GA execution
    for j = 1:number_execution
        fprintf('Executing GA for c = %d. Iteration = %d\n', c(idx), j);
        [nodes, solution, generations] = ga(G, c(idx), nNodes, population_size, q, num_elite, time_limit);
        results_ga{idx}(j, 1:c(idx)) = nodes;
        results_ga{idx}(j, end-1) = solution;
        results_ga{idx}(j, end) = generations;
        iterations_ga{idx}(j) = generations;

        % Check for new min or max
        if solution < local_min_ga
            local_min_ga = solution;
            min_nodes_ga{idx} = nodes;
        end
        if solution > local_max_ga
            local_max_ga = solution;
            max_nodes_ga{idx} = nodes;
        end
    end

    % Compute min, mean, max
    min_ga(idx) = local_min_ga;
    mean_ga(idx) = mean(results_ga{idx}(:, end-1));
    max_ga(idx) = local_max_ga;
end

% Prepare the header
fprintf('Method       | c value | Max value | Min value | Mean Value | Iterations/Generations | Nodes Min\n');
fprintf('-------------|---------|-----------|-----------|------------|------------------------|---------------------------------\n');

Display the results in a formatted way
for idx = 1:length(c)
    fprintf('GRASP        | %7d | %9.2f | %9.2f | %10.2f | %10.2f | %s\n', ...
        c(idx), max_grasp(idx), min_grasp(idx), mean_grasp(idx), mean(iterations_grasp{idx}), cellArrayToCSVString(min_nodes_grasp(idx)));
end

for idx = 1:length(c)
    fprintf('GA           | %7d | %9.2f | %9.2f | %10.2f | %10.2f | %s\n', ...
        c(idx), max_ga(idx), min_ga(idx), mean_ga(idx), mean(iterations_ga{idx}), cellArrayToCSVString(min_nodes_ga(idx)));
end

%--------------- GRASP ALGORITHM SECTION -------------------------

function [solution_nodes, solution, num_iter] = grasp(G, n, r, time_limit)
    % Initialize timer and counter
    t = tic;
    num_iter = 0;

    % Perform initial greedy randomized solution
    solution_nodes = greedyRandomized(G, n, r);
    solution = ConnectedNP(G, solution_nodes);

    % Main loop with time limit
    while toc(t) < time_limit
        % Generate new solution via the greedy randomized method
        current_nodes = solution_nodes;
        current_solution = solution;

        % Initialize local search parameters
        improved = true;
        while improved
            N = numnodes(G);
            aux = setdiff(1:N, current_nodes);
            best_neigh_solution = current_solution;
            best_curr_neigh = current_nodes;
            
            % Explore neighborhood
            for a = current_nodes
                for b = aux
                    neighbors = [setdiff(current_nodes,a), b];
                    neigh_solution = ConnectedNP(G, neighbors);
                    if neigh_solution < best_neigh_solution
                        best_neigh_solution = neigh_solution;
                        best_curr_neigh = neighbors;
                    end
                end
            end
            
            % Check if there is improvement
            if best_neigh_solution < current_solution
                current_solution = best_neigh_solution;
                current_nodes = best_curr_neigh;
            else
                improved = false;
            end
        end

        % Compare with global best solution
        if current_solution < solution
            solution = current_solution;
            solution_nodes = current_nodes;
        end

        num_iter = num_iter + 1; % Update iteration counter
    end
end

function s = greedyRandomized(G, n, r)
    
    E= 1:numnodes(G);
    s= [];
    for i= 1:n
        R= [];
        for j= E
            R= [R ; j ConnectedNP(G,[s j])];
        end
        R= sortrows(R,2);
        e= R(randi(r),1);
        s= [s e];
        E= setdiff(E,e);
    end

end

%--------------- GENETIC ALGORITHM SECTION ------------------------
function [solution_nodes, solution, num_generation] = ga(G, n, nNodes, population_size, q, num_elite, time_limit)
    num_generation = 0;
    %Generating the first population
    current_population = zeros([population_size, n+1]); 
    for i = 1:population_size
        current_individual = randperm(nNodes, n);
        current_population(i, 1:end-1) = current_individual;
    end

    num_generation = num_generation + 1;
    
    %initial fitness value computation
    for i = 1:population_size
        current_fitness = ConnectedNP(G, current_population(i, 1:end-1));
        current_population(i, end) = current_fitness;
    end
    
    current_population = sortrows(current_population, (n+1));
    solution = current_population(1, n+1);
    solution_nodes = current_population(1, 1:end-1);
    
    t = tic;

    while toc(t) < time_limit

        new_population = zeros(population_size, n+1);
        
        % Generate new individuals with crossover and mutation
        for i = 1:population_size
            new_individual = crossover_tournament(current_population, n);
            
            %veryfing mutation probability
            if rand < q
                new_individual = mutation(new_individual, nNodes, n);
            end
   
            new_fitness = ConnectedNP(G, new_individual(1:end-1)); %evaluating the fitness of the new indiviual
            new_individual(end) = new_fitness;
            new_population(i, :) = new_individual; %adding the new individual to the new population
        end
        
        % Sort the new population by fitness
        new_population = sortrows(new_population, (n+1));
        
        % Select the top m individuals as elite
        elite_individuals = new_population(1:num_elite, :);  
        
        % Combine current and new populations
        combined_population = [current_population; elite_individuals];
        combined_population = sortrows(combined_population, (n+1));
        
        % Keep the best num_of_individual from combined population
        current_population = combined_population(1:population_size, :);

        %i create a new generation
        num_generation = num_generation+1;

        % Track the best solution
        if current_population(1, end) < solution
            solution = current_population(1, end);
            solution_nodes = current_population(1, 1:end-1);
        end

    end
end

function new_individual = crossover_tournament(population, n)
    % Tournament selection to pick parent 1 (excluding the fitness
    % position)
    idx1 = randi(size(population - 1, 1));
    idx2 = randi(size(population - 1, 1));
    if population(idx1, n+1) > population(idx2, n+1)  %the last column is fitness
        parent_1 = population(idx1, 1:end-1);
    else
        parent_1 = population(idx2, 1:end-1);
    end

    % Tournament selection to pick parent 2
    idx1 = randi(size(population - 1, 1));
    idx2 = randi(size(population - 1, 1));
    if population(idx1, n+1) > population(idx2, n+1)
        parent_2 = population(idx1, 1:end-1);
    else
        parent_2 = population(idx2, 1:end-1);
    end

    aux = union(parent_1, parent_2);
    aux2 = randperm(length(aux), n);
    new_individual = aux(aux2);


    new_individual(n+1) = 0; % Temporary value, compute actual fitness later

    return
end

%mutation function
function mutated_individual = mutation(individual, nNodes, n)
   
    mutated_individual = individual;

    % selecting a random gene
    mutation_point = randi(n);

    % selecting a node that is not currently in the individual
    possible_nodes = setdiff(1:nNodes, individual(1:n));
    new_node = possible_nodes(randi(length(possible_nodes)));

    % applying mutation
    mutated_individual(mutation_point) = new_node;
    mutated_individual(n+1) = 0; % Imposta il fitness a 0 per indicare che deve essere ricalcolato

end





% ------------- Auxiliar function for dislaying the results ------------

% Function to convert cell arrays containing numeric arrays to comma-separated string
function nodeStr = cellArrayToCSVString(cellArray)
    numCells = numel(cellArray);
    nodeStr = '';
    for i = 1:numCells
        numericArray = cellArray{i};
        numericStr = ['[' num2str(numericArray) ']'];
        if i < numCells
            nodeStr = [nodeStr numericStr ', '];
        else
            nodeStr = [nodeStr numericStr];
        end
    end
end