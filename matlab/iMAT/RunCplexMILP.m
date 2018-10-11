function [milp_model] = RunCplexMILP (milp_model, print_level)
    % Creates and solves mixed-integer linear programming problems using the
    % Cplex format from the standard form of milp problem of a metabolic model

    % Initialize the CPLEX object
    try
        ILOGcplex = Cplex();
    catch ME
        error('CPLEX not installed or licence server not up')
    end

    % define problem
    ILOGcplex.Model.A = milp_model.S;
    ILOGcplex.Model.rhs = milp_model.rowub;
    ILOGcplex.Model.lhs = milp_model.rowlb;
    ILOGcplex.Model.ub = milp_model.ub;
    ILOGcplex.Model.lb = milp_model.lb;
    ILOGcplex.Model.obj = -milp_model.c;
    ind = find(-milp_model.c==1);
    vartype(milp_model.int_vars==1) = 'I';
    vartype(milp_model.int_vars==0) = 'C';
    % Make sure, that the vartype is in the correct orientation, cplex is quite picky here
    if size(vartype,1) > size(vartype,2)
        vartype = vartype';
    end
    ILOGcplex.Model.ctype = vartype;
    ILOGcplex.Start.x = zeros(size(milp_model.S,2),1);

    % parameters
    ILOGcplex.Param.mip.strategy.nodeselect.Cur = 0; % depth-first search
    ILOGcplex.Param.timelimit.Cur = 120;
    ILOGcplex.Param.output.clonelog.Cur = -1; % no log files

    % Optimize the problem
    ILOGcplex.solve();

    % get results
    milp_model.result = ILOGcplex.Solution;
    milp_model.result_vector = ILOGcplex.Solution.x;
    if ~isempty(ind)
        milp_model.result_opt = ILOGcplex.Solution.objval;
    else
        milp_model.result_opt = -ILOGcplex.Solution.objval;
    end
    milp_model.result_opt = round(milp_model.result_opt);

    if ILOGcplex.Solution.status==101
        milp_model.result_status = 0;
    else
        milp_model.result_status = 1;
    end
    milp_model.result_status_text = ILOGcplex.Solution.statusstring;

    if print_level > 0
        fprintf('\n*** RunCplexMILP ***\nOpt val: %d\nExit flag: %d\nExit text: %s\n', milp_model.result_opt, milp_model.result_status, milp_model.result_status_text);
    end
