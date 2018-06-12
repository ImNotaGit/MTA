function [milp_model] = RunGurobiMILP (milp_model, print_level)
    % Creates and solves mixed-integer linear programming problems using the
    % Gurobi format from the standard form of lp problem of a metabolic model

    % define problem
    neqind = find(milp_model.rowlb ~= milp_model.rowub);
    model.A = sparse([milp_model.S; milp_model.S(neqind,:)]);
    model.rhs = [reshape(milp_model.rowlb,length(milp_model.rowlb),1); reshape(milp_model.rowub(neqind),length(neqind),1)];
    model.sense(1:size(model.A,1)) = '=';
    model.sense(neqind) = '>';
    model.sense((size(milp_model.S,1)+1):size(model.A,1)) = '<';
    model.lb = milp_model.lb;
    model.ub = milp_model.ub;
    vartype(milp_model.int_vars==1) = 'I';
    vartype(milp_model.int_vars==0) = 'C';
    model.vtype = vartype;
    model.start = zeros(size(model.A,2),1);
    model.obj = double(-milp_model.c)+0;
    ind = find(-milp_model.c==1);
    model.modelsense = 'min';

    % parameter
    if print_level == 0
        params.OutputFlag = 0;
        params.DisplayInterval = 1;
    else
        params.OutputFlag = 1;
        params.DisplayInterval = 5;
    end
    params.Method = -1; % automatic
    params.TimeLimit =120;

    % call the solver
    resultgurobi = gurobi(model, params);

    % get results
    milp_model.result = resultgurobi;
    milp_model.result_vector = resultgurobi.x;
    if ~isempty(ind)
        milp_model.result_opt = resultgurobi.objval;
    else
        milp_model.result_opt = -resultgurobi.objval;
    end

    milp_model.result_status_text = resultgurobi.status;
    if strcmp(resultgurobi.status, 'OPTIMAL')
        milp_model.result_status = 0;
    else
        milp_model.result_status = 1;
    end

    if print_level > 0
        fprintf('\n*** RunGurobiMILP ***\nOpt val: %d\nExit text: %s\n', milp_model.result_opt, milp_model.result_status_text);
    end

end
