function [lp_model] = RunGurobiLP (lp_model, print_level)
    % Creates and solves linear programming problems using the
    % Gurobi format from the standard form of lp problem of a metabolic model

    % define problem
    neqind = find(lp_model.rowlb ~= lp_model.rowub);
    model.A = sparse([lp_model.S; lp_model.S(neqind,:)]);
    model.rhs = [reshape(lp_model.rowlb,length(lp_model.rowlb),1); reshape(lp_model.rowub(neqind),length(neqind),1)];
    model.sense(1:size(model.A,1)) = '=';
    model.sense(neqind) = '>';
    model.sense((size(lp_model.S,1)+1):size(model.A,1)) = '<';
    model.lb = lp_model.lb;
    model.ub = lp_model.ub;
    model.obj = double(-lp_model.c)+0;
    ind = find(-lp_model.c==1);
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
    params.Threads = 1; % single thread

    % call the solver
    resultgurobi = gurobi(model, params);

    % get results
    lp_model.result = resultgurobi;
    lp_model.result_vector = resultgurobi.x;
    if ~isempty(ind)
        lp_model.result_opt = resultgurobi.objval;
    else
        lp_model.result_opt = -resultgurobi.objval;
    end

    lp_model.result_status_text = resultgurobi.status;
    if strcmp(resultgurobi.status, 'OPTIMAL')
        lp_model.result_status = 0;
    else
        lp_model.result_status = 1;
    end

    if print_level > 0
        fprintf('\n*** RunGurobiLP ***\nOpt val: %d\nExit text: %s\n', lp_model.result_opt, lp_model.result_status_text);
    end

end
