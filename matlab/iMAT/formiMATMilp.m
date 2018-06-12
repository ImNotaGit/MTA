function [milp,v_act_fwr,v_act_bck,v_inact] = formiMATMilp(model, exp_data, solver)

% Create a MILP problem from the model and the gene expression/reaction data
    
    Private_DefineParameters;

    % INIT REACTION DATA
    reaction_data = zeros(length(model.lb),1);

    % Initialize reaction data (if not given) from gene expression data.
    if ~isempty(exp_data)
        for i = 1:length(model.rules)
            reaction_data(i) = evalExpRule2(model.rules{i}, exp_data);
        end
    end 
    % remove integers for dead end rxns
    tmp_reaction_data = reaction_data;
    for (i=1:length(reaction_data))
        if ((model.lb(i) == 0) && (model.ub(i) == 0) && (reaction_data(i) ~= 0))
            reaction_data(i) = 0;
        end
    end
    
    % Printout - Up/Down/Non expressed genes and reactions
    fprintf('Model gene expression: up=%d, Down=%d, Other=%d\n', sum(exp_data==1), sum(exp_data==-1), sum(exp_data==0));
    fprintf('Model reaction expression: up=%d, Down=%d, Other=%d\n', sum(reaction_data==1), sum(reaction_data==-1), sum(reaction_data==0));

    
    % CREATE MILP PROBLEM
    milp = model;
    milp.rowlb = model.b;   % rowlb/ub are the constraints on S*v. Initialize with the original constraints (b is 0 - steady state),
    milp.rowub = model.b;   % and then expand for the MILP constraints on the new variables - see below.    
    

    % ACTIVE REACTIONS - FORWARD
    % Expand S for the forward direction of active reactions, we want Vi to
    % be greater than epsilon. Note: All reactions have forward direction.
    % The constraint is Vi + Y(+)i *  (Vmin - epsilon) >= Vmin
    % so that: if Y(+)i == 0, Vi >= Vmin, if Y(+)i == 1, Vi >= epsilon
    v = find(reaction_data==1); n = length(v);  % v - indices of all active reactions (uni-directional and bi-directional)
    v_act_fwr = v;
    % m1 = rows: active reactions, set 1 for each reaction in the appropriate slot (the reaction index)
    % i.e. if the active reactions are 2,5,8 and there are 10 reactions, the matrix will look like
    % 0100000000
    % 0000100000
    % 0000000100
    m1 = sparse( [1:n], v, ones(n,1), n, size(milp.S,2) );  % size(milp.S,2) = all rxns
    m2 = eye(n) * (-FLUX_BOUND - ACTIVE_FLUX);
    milp.S = [milp.S, sparse(size(milp.S,1), n) ; m1 m2];   % size(milp.S,1) = all mets
    milp.rowlb = [milp.rowlb ; zeros(n,1)-FLUX_BOUND];
    milp.rowub = [milp.rowub ; zeros(n,1)+FLUX_BOUND];
    milp.var_ind = v;

    % ACTIVE REACTIONS - BACKWARD
    % Expand S for the backward direction of active reactions, we want Vi to
    % be lower than -epsilon. Note: Only bi-directional reactions have
    % backward direction.
    v = find(reaction_data==1 & model.lb<0);   n = length(v);
    v_act_bck = v;
    m1 = sparse( [1:n], v, ones(n,1), n, size(milp.S,2) );  
    m2 = eye(n) * (FLUX_BOUND + ACTIVE_FLUX);
    milp.S = [milp.S, sparse(size(milp.S,1), n) ; m1 m2];
    milp.rowlb = [milp.rowlb ; zeros(n,1)-FLUX_BOUND];
    milp.rowub = [milp.rowub ; zeros(n,1)+FLUX_BOUND];
    milp.var_ind = [milp.var_ind; v];

    % INACTIVE REACTIONS - FORWARD/BACKWARD
    % Expand S for inactive reactions, we want Vi to be 0
    % The constraint is Vmin * (1-Yi) <= Vi <= Vmax * (1-Yi)
    % so that: if Yi == 0, Vmax >= Vi >= Vmin, if Yi == 1, Vi = 0
    % Note that in practice we allow epsilon (currently
    % 0.1), so that -epsilon <= Vi <= epsilon. Different epsilons are used
    % for active and inactive reactions.
    [met_num, rxn_num] = size(model.S);
    v = find(reaction_data==-1); n = length(v);
    v_inact = v;
    m1 = sparse( [1:n], v, ones(n,1), n, size(milp.S,2) );
    m2 = eye(n) * (FLUX_BOUND - INACTIVE_FLUX );

    v2 = find(reaction_data==-1 & model.lb<0); n2 = length(v2);
    m3 = sparse( [1:n2], v2, ones(n2,1), n2, size(milp.S,2) );
    m4 = sparse( [1:n2], find(model.lb(v)<0), 1, n2, n) * (-FLUX_BOUND + INACTIVE_FLUX);

    milp.S = [milp.S, sparse(size(milp.S,1), n) ; m1 m2 ; m3 m4];
    milp.rowlb = [milp.rowlb ; zeros(n,1)-FLUX_BOUND ; zeros(n2,1)-FLUX_BOUND ];
    milp.rowub = [milp.rowub ; zeros(n,1)+FLUX_BOUND ; zeros(n2,1)+FLUX_BOUND];
    milp.var_ind = [milp.var_ind; v];
    milp.var_ind = [milp.var_ind; v2];
    
    % 
    n = size(milp.S,2) - rxn_num;
    milp.c = [zeros(rxn_num,1); ones(n, 1)];
    milp.int_vars =  milp.c;
    milp.lb = [milp.lb ; zeros(n, 1)];
    milp.ub = [milp.ub ; ones(n, 1)];
    milp.rxn_num = rxn_num;
    milp.met_num = met_num; 
    milp.gene_num = length(model.genes);
    milp.gene_exp = exp_data;
    milp.reaction_exp = tmp_reaction_data;
    milp.new_reaction_exp = reaction_data;
    
    % find opt solution value
    switch solver
        case 'tomlab'
            wt_result = RunTomlabMILP(milp,1);
        case 'cplex'
            wt_result = RunCplexMILP(milp,1);
        case 'gurobi'
            wt_result = RunGurobiMILP(milp,1);
    end
    if  (wt_result.result_status ~= 0)
        error('ERROR: failed running wt optimization');
    end
    % save wt result
    wt_flux = wt_result.result_vector(1:rxn_num);
    milp.wt_flux = wt_flux;
    milp.start_sol = wt_result.result_vector;
    milp.wt_opt = round(wt_result.result_opt);
    milp.reaction_data = reaction_data;
    
end
