function[success, unsuccess, wrong, dvs, dvu, dvw, diff] = calculateDiff(v_ref,fwd,bck,cons_rxns_fb,result)

fwd = intersect(fwd,cons_rxns_fb);
bck = intersect(bck,cons_rxns_fb);
both = intersect(fwd,bck);

[u,ia,ib] = intersect(fwd,both);
fwd(ia) = [];
[u,ia,ib] = intersect(bck,both);
bck(ia) = [];

fwd_wrong = fwd(find(v_ref(fwd)<0 & result(fwd)>abs(v_ref(fwd))));
bck_wrong = bck(find(v_ref(bck)>0 & result(bck)<-(v_ref(bck))));
wrong = [fwd_wrong; bck_wrong];
%diff_fwd_wrong = sum(abs(abs(v_ref(fwd_wrong)) - abs(result(fwd_wrong))));
%diff_bck_wrong = sum(abs(abs(v_ref(bck_wrong)) - abs(result(bck_wrong))));
dvw = abs(abs(v_ref(wrong)) - abs(result(wrong)));
diffw = sum(dvw);

[u,ia,ib] = intersect(fwd,fwd_wrong);
fwd(ia) = [];
[u,ia,ib] = intersect(bck,bck_wrong);
bck(ia) = [];

fwd_success = fwd(find(result(fwd)>v_ref(fwd)));
bck_success = bck(find(result(bck)<v_ref(bck)));
success = [fwd_success; bck_success; both];
%diff_fwd_success = sum(abs(v_ref(fwd_success)-result(fwd_success)));
%diff_bck_success = sum(abs(v_ref(bck_success)-result(bck_success)));
%diff_both = sum(abs(v_ref(both)-result(both)));
dvs = abs(v_ref(success)-result(success));
diffs = sum(dvs);

fwd_unsuccess = fwd(find(result(fwd)<v_ref(fwd)));
bck_unsuccess = bck(find(result(bck)>v_ref(bck)));
unsuccess = [fwd_unsuccess; bck_unsuccess];
%diff_fwd_unsuccess = sum(abs(v_ref(fwd_unsuccess)-result(fwd_unsuccess)));
%diff_bck_unsuccess = sum(abs(v_ref(bck_unsuccess)-result(bck_unsuccess)));
dvu = abs(v_ref(unsuccess)-result(unsuccess));
diffu = sum(dvu);

%diff = diff_fwd_sucess + diff_bck_sucess + diff_both - diff_fwd_unsucess - diff_bck_unsucess - diff_fwd_wrong - diff_bck_wrong;
diff = diffs - diffu - diffw;