function[diff] = calculateDiff(v_ref,fwd,bck,cons_rxns_fb,result)

fwd = intersect(fwd,cons_rxns_fb);
bck = intersect(bck,cons_rxns_fb);
both = intersect(fwd,bck);

[u,ia,ib] = intersect(fwd,both);
fwd(ia) = [];
[u,ia,ib] = intersect(bck,both);
bck(ia) = [];

fwd_wrong = fwd(find(v_ref(fwd)<0 & result(fwd)>abs(v_ref(fwd))));
bck_wrong = bck(find(v_ref(bck)>0 & result(bck)<-(v_ref(bck))));

diff_fwd_wrong = sum(abs(abs(v_ref(fwd_wrong)) - abs(result(fwd_wrong))));
diff_bck_wrong = sum(abs(abs(v_ref(bck_wrong)) - abs(result(bck_wrong))));

[u,ia,ib] = intersect(fwd,fwd_wrong);
fwd(ia) = [];
[u,ia,ib] = intersect(bck,bck_wrong);
bck(ia) = [];

fwd_sucess = fwd(find(result(fwd)>v_ref(fwd)));
fwd_unsucess = fwd(find(result(fwd)<v_ref(fwd)));
bck_sucess = bck(find(result(bck)<v_ref(bck)));
bck_unsucess = bck(find(result(bck)>v_ref(bck)));

diff_fwd_sucess = sum(abs(v_ref(fwd_sucess)-result(fwd_sucess)));
diff_fwd_unsucess = sum(abs(v_ref(fwd_unsucess)-result(fwd_unsucess)));
diff_bck_sucess = sum(abs(v_ref(bck_sucess)-result(bck_sucess)));
diff_bck_unsucess = sum(abs(v_ref(bck_unsucess)-result(bck_unsucess)));
diff_both = sum(abs(v_ref(both)-result(both)));

diff = diff_fwd_sucess + diff_bck_sucess + diff_both - diff_fwd_unsucess - diff_bck_unsucess - diff_fwd_wrong - diff_bck_wrong;