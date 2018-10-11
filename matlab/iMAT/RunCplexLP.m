function [lp_model] = RunCplexLP (lp_model, print_level)
    % Creates and solves linear programming problems using the
    % Cplex format from the standard form of lp problem of a metabolic model
    
    % Initialize the CPLEX object
    try
        ILOGcplex = Cplex();
    catch ME
        error('CPLEX not installed or licence server not up')
    end

    % define problem
    ILOGcplex.Model.sense = 'minimize';
    ILOGcplex.Model.obj   = -lp_model.c;
    ILOGcplex.Model.lb    = lp_model.lb;
    ILOGcplex.Model.ub    = lp_model.ub;
    ILOGcplex.Model.A     = lp_model.S;
    ILOGcplex.Model.lhs   = lp_model.rowlb;
    ILOGcplex.Model.rhs   = lp_model.rowub;
    ind = find(-lp_model.c==1);
    
    % parameters
    if print_level==0
        ILOGcplex.DisplayFunc=[];
    end
    ILOGcplex.Param.threads.Cur = 1; % single thread
    ILOGcplex.Param.output.clonelog.Cur = -1; % no log files

    % Optimize the problem
    ILOGcplex.solve();

    % get results
    lp_model.result = ILOGcplex.Solution;
    lp_model.result_vector = ILOGcplex.Solution.x;
    if ~isempty(ind)
        lp_model.result_opt = ILOGcplex.Solution.objval;
    else
        lp_model.result_opt = -ILOGcplex.Solution.objval;
    end

    if ILOGcplex.Solution.status==1
        lp_model.result_status = 0;
    else
        lp_model.result_status = 1;
    end
    lp_model.result_status_text = ILOGcplex.Solution.statusstring;

    if print_level > 0
        fprintf('\n*** RunCplexLP ***\nOpt val: %d\nExit flag: %d\nExit text: %s\n', lp_model.result_opt, lp_model.result_status, lp_model.result_status_text);
    end

end
