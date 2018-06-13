function warmup_points = model_sampling_warmup(model, lp_solver, points_count, is_verbose)
    if (nargin < 3)
        points_count = 5000;
    end
    if (nargin < 4)
        is_verbose = false;
    end

    rxn_count = length(model.rxns);
    
    % Create orthogonal warmup points
    warmup_points_orth = get_orth_points(model, is_verbose, lp_solver);
    orth_points_count = length(warmup_points_orth);
    
    warmup_points = zeros(rxn_count, points_count);

    orth_point_order = randperm(orth_points_count);

    % Combine point sets
    if is_verbose
        fprintf('start generating warmup points\n')
    end
    parfor i = 1:points_count
        if (i <= orth_points_count)
            % Ensure that each direction is used at least once
            orth_point_idx = orth_point_order(i);
        else
            % All direction already used
            orth_point_idx = ceil(rand * orth_points_count);
        end
        orth_point = warmup_points_orth(:, orth_point_idx);
        random_point = get_random_point(model, lp_solver);
        r = rand;
        the_point = orth_point * r + random_point * (1 - r);
        warmup_points(:, i) = the_point;
        % if paralleling, print_progress does not print in real-time, so disabled
        %if is_verbose && mod(i, 100) == 0
        %    print_progress(i / points_count);
        %end
    end
    if is_verbose
        fprintf('done generating warmup points\n')
    end
end

function point = get_random_point(model, lp_solver)
    % Create random objective function
    c = rand(length(model.rxns), 1) - 0.5;
    point = get_opt_point(model, c, lp_solver);
end

function points = get_orth_points(model, is_verbose, lp_solver)
    if is_verbose
        fprintf('start generating orthogonal warmup points\n');
    end
    rxns_count = length(model.rxns);
    points1 = zeros(rxns_count, rxns_count);
    points2 = zeros(rxns_count, rxns_count);
    parfor rxn_idx = 1:rxns_count
        % Pick the next flux to optimize, cycles though each reaction
        % alternates minimization and maximization for each cycle
        c = zeros(size(model.c));
        c(rxn_idx) = -1;
        point = get_opt_point(model, c, lp_solver);
        points1(:, rxn_idx) = point;
        c(rxn_idx) = 1;
        point = get_opt_point(model, c, lp_solver);
        points2(:, rxn_idx) = point;
        % if paralleling, print_progress does not print in real-time, so disabled
        %if is_verbose && mod(rxn_idx, 100) == 0
        %    print_progress(rxn_idx / rxns_count);
        %end
    end
    points = [points1 points2];
    fprintf('done generating orthogonal warmup points\n');
end

function point = get_opt_point(model, c, lp_solver)
    c = c / norm(c);
    model.c = c;
    % Find optimal solution
    switch lp_solver
        case 'tomlab'
            result = RunTomlabLP(model, 0);
        case 'cplex'
            result = RunCplexLP(model, 0);
        case 'gurobi'
            result = RunGurobiLP(model, 0);
    end
    point = result.result_vector;
end
