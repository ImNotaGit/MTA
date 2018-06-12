function [milp_model] = RunTomlabMILP (milp_model,print)
% Creates and solves mixed-integer linear programming problems using the
% TOMLAB format from the standard form of milp problem of a metabolic model
Name=[];%milp_model.model_name;
% Problem formulated as a maximum problem
A = milp_model.S;
b_U = milp_model.rowub;
b_L = milp_model.rowlb;
c = milp_model.c;
c=-c;
ind = find(c==1);
ind2 = find(c==-1);
[m,n] = size(A);
x_L = milp_model.lb;
x_U = milp_model.ub;
x_0 = zeros(n,1);
%fprintf('MILP problem. Variables %d. Knapsacks %d\n',n,m);
IntVars = milp_model.int_vars;
x_min = x_L; x_max = x_U; f_Low = -1E7; % f_Low <= f_optimal must hold
f_opt = [];%Optimal function value(s), if known (Stationary points) not used
nProblem = []; % Problem number not used
fIP = []; % Do not use any prior knowledge
xIP = []; % Do not use any prior knowledge
setupFile = []; % Just define the Prob structure, not any permanent setup file
x_opt = []; % The optimal integer solution is not known
VarWeight = []; % No variable priorities, largest fractional part will be used
KNAPSACK = 0; % Run without the knapsack heuristic
% Assign routine for defining a MIP problem.
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, ...
    nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.MIP.cpxControl.TILIM = 120;
Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
Prob.Solver.Alg = 2; % Depth First, then Breadth search
% Calling driver routine tomRun to run the solver.
% The 1 sets the print level after optimization.
milp_model.tomlab_result = tomRun('cplex', Prob, print);
milp_model.result_vector = milp_model.tomlab_result.x_k;
if ~isempty(ind)
    lp_model.result_opt = milp_model.tomlab_result.f_k;
else
    lp_model.result_opt = -milp_model.tomlab_result.f_k;
end
milp_model.result_opt = round(-milp_model.tomlab_result.f_k);
milp_model.result_status = milp_model.tomlab_result.ExitFlag;
milp_model.result_status_text = milp_model.tomlab_result.ExitText;
%fprintf('Exit flag: %d\nExit text: %s\n',milp_model.result_status,milp_model.result_status_text);
end