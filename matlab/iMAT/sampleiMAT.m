function[sample_points,v_ref] = sampleiMAT(model, discrete_data, milp_solver, lp_solver)
% function for running iMAT
% model: the metabolic model
% discrete_data: vector of discretized gene expression of the source state, i.e. the output from get_binary_expH.m

Private_DefineParameters;

[milp,v_act_fwr,v_act_bck,v_inact] = formiMATMilp(model, discrete_data, milp_solver);

rxns = [v_act_fwr; v_act_bck; v_inact];
integers = milp.start_sol(length(model.rxns)+1:end);

int_act_fwr = integers(1:length(v_act_fwr));
int_act_bck = integers(length(v_act_fwr)+1:length(v_act_fwr)+length(v_act_bck));
int_inact = integers(length(v_act_fwr)+length(v_act_bck)+1:end);

active_fwr  = v_act_fwr(find(int_act_fwr==1));
active_bck  = v_act_bck(find(int_act_bck==1));
inactive  = v_inact(find(int_inact==1));
inactive_rev = find(model.lb(inactive)<0);
inactive_non_rev = find(model.lb(inactive)==0);

model.lb(active_fwr) = ACTIVE_FLUX;
model.ub(active_bck) = -ACTIVE_FLUX;
model.lb(inactive(inactive_rev)) = -INACTIVE_FLUX;
model.lb(inactive(inactive_non_rev)) = 0;
model.ub(inactive) = INACTIVE_FLUX;

[sample_points]= model_sample(model, 2000, lp_solver); % 2000 is the times of repeated iMAT runs (sampling) for each source state sample
v_ref = mean(sample_points,2); % mean of samples