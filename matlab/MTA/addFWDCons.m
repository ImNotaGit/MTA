%Construct the FWD constraints
function[model] = addFWDCons(rxns,ref_val,model,EPS)

DEFINE_PARAM;
[m,n] = size(model.S);
col1 = zeros(m,1);
col2 = zeros(m,1);
row1 = zeros(1,n+2);
row2 = zeros(1,n+2);
row1([rxns,n+1,n+2]) = [1, -(ref_val+EPS),-V_MIN]; %v_i - y_i_F(v_i_ref+eps) - y_i*v_i_min
row2([n+1,n+2]) = [1, 1]; %y_i_F + y_i
model.S = [model.S col1 col2];
model.S = [model.S;row1;row2];
model.rowlb = [model.rowlb; 0; 1];
model.rowub = [model.rowub; V_MAX_C; 1];
model.lb = [model.lb;0;0];
model.ub = [model.ub;1;1];