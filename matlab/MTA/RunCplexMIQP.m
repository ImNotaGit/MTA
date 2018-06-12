function [qp_model] = RunCplexMIQP (qp_model, print_level)
    %Running a MIQP problem according to Cplex formulation

    % Initialize the CPLEX object
    try
        ILOGcplex = Cplex();
    catch ME
        error('CPLEX not installed or licence server not up')
    end

    % define problem
    ILOGcplex.Model.A = qp_model.S;
    ILOGcplex.Model.rhs = qp_model.rowub;
    ILOGcplex.Model.lhs = qp_model.rowlb;
    ILOGcplex.Model.ub = qp_model.ub;
    ILOGcplex.Model.lb = qp_model.lb;
    ILOGcplex.Model.obj = qp_model.c;
    ILOGcplex.Model.Q = qp_model.F;
    vartype(qp_model.int_vars==1) = 'I';
    vartype(qp_model.int_vars==0) = 'C';
    % Make sure, that the vartype is in the correct orientation, cplex is quite picky here
    if size(vartype,1) > size(vartype,2)
        vartype = vartype';
    end
    ILOGcplex.Model.ctype = vartype;
    ILOGcplex.Start.x = zeros(size(qp_model.S,2),1);

    % parameters
    ILOGcplex.Param.timelimit.Cur = 120;

    % Optimize the problem
    ILOGcplex.solve();

    % get results
    qp_model.result = ILOGcplex.Solution;
    qp_model.result_vector = ILOGcplex.Solution.x;
    qp_model.result_opt = ILOGcplex.Solution.objval;

    qp_model.result_status = ILOGcplex.Solution.status;
    qp_model.result_status_text = ILOGcplex.Solution.statusstring;

    if print_level > 0
        fprintf('\n*** RunCplexMIQP ***\nOpt val: %d\nExit flag: %d\nExit text: %s\n', qp_model.result_opt, qp_model.result_status, qp_model.result_status_text);
    end
