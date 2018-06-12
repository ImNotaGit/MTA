%Construct the BCK constraints
function[model] = addBCKCons(rxns,ref_val,model,EPS)

DEFINE_PARAM;
[m,n] = size(model.S);
col1 = zeros(m,1);
col2 = zeros(m,1);
row1 = zeros(1,n+2);
row2 = zeros(1,n+2);
row1([rxns,n+1,n+2]) = [1, -(ref_val-EPS),-V_MAX]; %v_i - y_i_B(v_i_ref-eps) - y_i*v_i_max
row2([n+1,n+2]) = [1, 1]; %y_i_B + y_i
model.S = [model.S col1 col2];
model.S = [model.S;row1;row2];
model.rowlb = [model.rowlb; V_MIN_C; 1];
model.rowub = [model.rowub; 0; 1];
model.lb = [model.lb;0;0];
model.ub = [model.ub;1;1];
