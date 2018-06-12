function[milp,v_act_fwr,v_act_bck,v_inact] = testmymilp(model, discrete_data, milp_solver)
% function for running iMAT
% model: the metabolic model
% discrete_data: vector of discretized gene expression of the source state, i.e. the output from 01_get_binary_expH.m

Private_DefineParameters;

[milp,v_act_fwr,v_act_bck,v_inact] = formiMATMilp(model, discrete_data, milp_solver);
