function [discrete_rxns_vector, geneState]=createDiscreteRxns_human(source, target, smean, tmean, model, stat)
% function for determining the directions of reaction fluxes changes from the "source state" to the "target state", producing a vector of -1/0/1's (down/no change/up)
% source: vector/matrix of gene expression in the source state
% target: vector/matrix of gene expression in the target state
% model: the metabolic model
% stat==0: for a pair of sample; stat==1: for two groups of samples, each group with multiple samples, in which case source and target are matrices where the columns are samples, and smean and tmean are the mean for source and target samples, respectively

prc = 0.96 % for stat==0, try different values to get about 100 reactions with -1/+1 direction in the output, for the efficiency of MTA
pval = 0.007 % for stat==1, try different values to get about 100 reactions with -1/+1 direction in the output, for the efficiency of MTA
geneState = zeros(size(model.genes_unique));

if stat ==1
    ignoreG = zeros(size(model.genes_unique));

    for i = 1:length(model.genes_unique)
        if length(find(model.rxnGeneMat(:,i)))>10 % if too many reactions mapped to a gene, then ignore this gene # always check manually!!!! and manually check everything before running. But, apparently this ignoreG is calculated but never actually used!
            ignoreG(i)=1;
        end
    end

    for i = 1:length(model.genes_unique)
        [h1(i),p1(i)] = ttest(target(i,:),smean(i),0.05,'left'); % ttest2(target(i,:),source(i,;),0.05, 'left'); for new version of matlab
    end

    for i = 1:length(model.genes_unique)
        [h2(i),p2(i)] = ttest(target(i,:),smean(i),0.05,'right');
    end
    
    discrete_rxns_vector = zeros(size(model.rxns));

    for i = 1:length(model.genes_unique)
        if p2(i)<=pval
            geneState(i) = 1;
        end
        if p1(i)<=pval
            geneState(i) = -1;
        end
    end
end

if stat ==0
    diff = abs(source - target);
    thr = quantile(diff , prc);
    
    geneState(source - target>thr) = -1;
    geneState(target - source>thr) = 1;
    
end

geneState2 = zeros(size(model.genes));

for i=1:length(model.genes_unique)
    ind = model.genes_unique_map==i;
    geneState2(ind) = geneState(i);
end
    
for i = 1:length(model.rules)
    discrete_rxns_vector(i) = evalExpRule(model.rules{i}, geneState2);
end
