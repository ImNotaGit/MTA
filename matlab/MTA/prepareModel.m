function[model] = prepareModel

load ecoli_model
model = ecoli_model;
model.rowlb = zeros(length(model.mets),1);
model.rowub = zeros(length(model.mets),1);
model.int_vars = zeros(length(model.rxns),1);

model.lb(model.lb==-1000) = -50;
model.ub(model.ub==1000) = 50;

ex_rxns = strmatch('EX_',model.rxns);
media=find(model.lb(ex_rxns)<0);
model.lb(ex_rxns(media)) = -5;
Res = RunTomlabLP(model,1);
model.lb(find(model.c==1)) = 0.2*Res.result_opt;