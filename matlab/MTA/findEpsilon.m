function[p1,p2] = findEpsilon(sample_points_source, v_ref, discrete_rxns_vector,eps)

e_rxns = find(discrete_rxns_vector==1 | discrete_rxns_vector==2);
r_rxns = find(discrete_rxns_vector==-1 | discrete_rxns_vector==2);

for i=1:length(e_rxns)
    [h1(i),p1(i)] = ttest(sample_points_source(e_rxns(i),:),v_ref(e_rxns(i))+eps,0.05,'left');
end

for i=1:length(r_rxns)
    [h2(i),p2(i)] = ttest(sample_points_source(r_rxns(i),:),v_ref(r_rxns(i))-eps,0.05,'right');
end

(length(find(p1<=0.05))+length(find(p2<=0.05)))/(length(p1)+length(p2))