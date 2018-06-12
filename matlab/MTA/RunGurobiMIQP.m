function [qp_model] = RunGurobiMIQP (qp_model, print_level)
    %Running a MIQP problem according to Gurobi formulation

    % define problem
    neqind = find(qp_model.rowlb ~= qp_model.rowub);
    model.A = sparse([qp_model.S; qp_model.S(neqind,:)]);
    model.rhs = [reshape(qp_model.rowlb,length(qp_model.rowlb),1); reshape(qp_model.rowub(neqind),length(neqind),1)];
    model.sense(1:size(model.A,1)) = '=';
    model.sense(neqind) = '>';
    model.sense((size(qp_model.S,1)+1):size(model.A,1)) = '<';
    model.lb = qp_model.lb;
    model.ub = qp_model.ub;
    vartype(qp_model.int_vars==1) = 'I';
    vartype(qp_model.int_vars==0) = 'C';
    model.vtype = vartype;
    model.start = zeros(size(model.A,2),1);
    model.obj = double(qp_model.c)+0;
    model.modelsense = 'min';
    model.Q = sparse(0.5*qp_model.F); % tomlab/cplex optimizes 0.5x'Fx + c'x, while gurobi optimizes x'Qx + c'x

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
    params.Threads = 1; % single thread

    % call the solver
    resultgurobi = gurobi(model, params);

    % get results
    qp_model.result = resultgurobi;
    qp_model.result_vector = resultgurobi.x;
    qp_model.result_opt = resultgurobi.objval;

    qp_model.result_status_text = resultgurobi.status;
    if strcmp(resultgurobi.status, 'OPTIMAL')
        qp_model.result_status = 0;
    else
        qp_model.result_status = 1;
    end

    if print_level > 0
        fprintf('\n*** RunGurobiMIQP ***\nOpt val: %d\nExit text: %s\n', qp_model.result_opt, qp_model.result_status_text);
    end

end
