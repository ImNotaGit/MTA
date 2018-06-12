%Running a MIQP problem according to Tomlab formulation
function [qp_model] = RunTomlabMIQP (qp_model, print_level)

Name = qp_model.description;
A = qp_model.S;
b_U = qp_model.rowub;
b_L = qp_model.rowlb;
c= qp_model.c;
[m,n] = size(A);
x_L = qp_model.lb;
x_U = qp_model.ub;
x_0 = zeros(n,1);
IntVars = qp_model.int_vars;
if (print_level > 0)
    fprintf('qp problem. Variables %d. Knapsacks %d\n',n,m);
end
x_min = x_L; x_max = x_U;
F = qp_model.F;
Prob = miqpAssign(F, c, A, b_L, b_U, x_L, x_U, [],IntVars,[],[],[], Name,[],[]);
Prob.MIP.cpxControl.TILIM = 120;
qp_model.tomlab_result = tomRun('cplex', Prob, print_level);
qp_model.result_vector = qp_model.tomlab_result.x_k;
qp_model.result_opt = qp_model.tomlab_result.f_k;
qp_model.result_status = qp_model.tomlab_result.ExitFlag;
qp_model.result_status_text = qp_model.tomlab_result.ExitText;
if print_level > 0
    fprintf('\n*** RunTomlabQP ***\nOpt val: %d\nExit flag: %d\nExit text: %s\n',qp_model.result_opt,qp_model.result_status,qp_model.result_status_text);
end

