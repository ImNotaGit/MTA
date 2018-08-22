function [bin] = get_binary_expH(model, exp, stat)
% function for discretizing gene expressions of the "source state" into low/medium/high (-1/0/1) bins
% model: the metabolic model
% exp: vector/matrix of expression levels of the genes in the source state, with the same length and in the same order as model.genes_unique_names, 1496 for recon1 (model.mat)
% stat==0: discretize for one sample based on 25% and 75% quantile; stat==1: Noam mentioned that it's for multiple samples where exp is a matrix, but after examining the code it seems that it's just another way of discretizing for one sample: based on mean +- 0.7sd for the greater than 0 values.

if stat==0
    top = quantile(exp,0.75);
    down = quantile(exp, 0.25);

    tmp = zeros(size(model.genes_unique));
    tmp(exp>top) = 1;
    tmp(exp<down) = -1;

    bin = zeros(size(model.genes));
end

if stat==1
    model_exp = exp;

    meanVal = mean(model_exp(model_exp>0))
    dev = std(model_exp(model_exp>0))
    tmp = zeros(size(model.genes));
    0.5*dev

    for i=1:length(model_exp)
        if (model_exp(i))>meanVal+0.7*dev
            tmp(i)=1;
        end
        if (model_exp(i))<meanVal-0.7*dev
            tmp(i)=-1;
        end
    end
end

for i=1:length(model.genes_unique)
    ind = model.genes_unique_map==i;
    bin(ind) = tmp(i);
end

