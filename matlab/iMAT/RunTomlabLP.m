function [lp_model] = RunTomlabLP (lp_model, print_level)
    % Creates and solves linear programming problems using the
    % TOMLAB format from the standard form of lp problem of a metabolic model
    %Name=lp_model.name;
    %Name =lp_model.model_name;
    Name = lp_model.description;
    % Problem formulated as a maximum problem
    A = lp_model.S;
    b_U = lp_model.rowub;
    b_L = lp_model.rowlb;
    c = lp_model.c;
    c=-c;
    ind = find(c==1);
    ind2 = find(c==-1);
    
    [m,n] = size(A);
    x_L = lp_model.lb;
    x_U = lp_model.ub;
    x_0 = zeros(n,1);
    if (print_level > 0)
        fprintf('lp problem. Variables %d. Knapsacks %d\n',n,m);
    end
    x_min = x_L; x_max = x_U; 
    
    % Formulating the lp problem struct
    Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, [], [], [], x_min, x_max, [], []);

    % Calling driver routine tomRun to run the solver.
    lp_model.tomlab_result = tomRun('cplex', Prob, print_level);
    lp_model.result_vector = lp_model.tomlab_result.x_k;
    if ~isempty(ind)
    lp_model.result_opt = lp_model.tomlab_result.f_k;
    else
        lp_model.result_opt = -lp_model.tomlab_result.f_k;
    end
        
    lp_model.result_status = lp_model.tomlab_result.ExitFlag;
    lp_model.result_status_text = lp_model.tomlab_result.ExitText;
    if print_level > 0
        fprintf('\n*** RunTomlabLP ***\nOpt val: %d\nExit flag: %d\nExit text: %s\n',lp_model.result_opt,lp_model.result_status,lp_model.result_status_text);
    end
end
