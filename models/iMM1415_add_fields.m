load iMM1415
load iMM1415_rules % rules, created in R, see the beginning of ../Rwrapper/utils.R

 % the empty elements in rules is displayed as [1x0 char], which is equivalent to '' by strcmp, but I explicitly make them ''
rules(find(strcmp(rules, ''))) = {''};
rules
iMM1415.rules = rules;
iMM1415.rowlb = iMM1415.b; % all 0
iMM1415.rowub = iMM1415.b; % all 0
iMM1415
save iMM1415 iMM1415
