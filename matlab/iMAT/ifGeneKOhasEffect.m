function [hasEffect, constrRxnNames] = ifGeneKOhasEffect(model, geneList, downRegFraction)
% Deletes one gene and constrain the reactions affected to zero and appends '_deleted' to the gene(s)
%
% INPUT:
%    model:             COBRA model with the appropriate constrains for a
%                       particular condition
%    geneList:          A single gene to be deleted
%    downRegFraction:   Fraction of the original bounds that the reactions
%                       corresponding to downregulated genes will be assigned
%                       (Default = 0 corresponding to a full deletion)
%
% OUTPUTS:
%    model:             COBRA model with the selected genes deleted
%    hasEffect:         True if the gene deletion has an effect on the model
%    constrRxnNames:    Reactions that are associated to the genes in `geneList`

if (nargin < 3)
    downRegFraction = 0;
end

if (~iscell(geneList))
    geneName = geneList;
    clear geneList;
    geneList{1} = geneName;
end

if (~isfield(model,'genes'))
    error('Gene-reaction associations not included with the model');
end

%RxnGeneMat is required for this function, so we will have to build it if
%it does not exist
if ~isfield(model,'rxnGeneMat')
    model = buildRxnGeneMat(model);
end

hasEffect = false;
constrRxnNames = {};

% Find gene indices in model
[isInModel,geneInd] = ismember(geneList, model.genes);

if (isInModel)
  
  %mark genes for deletion
  model.genes(geneInd) = strcat(model.genes(geneInd),'_deleted');

    % Find rxns associated with this gene
    rxnInd = find(any(model.rxnGeneMat(:,geneInd),2));
    if (~isempty(rxnInd))
        x = true(size(model.genes));
        % set genes marked "_deleted" to false
        x(~cellfun('isempty',(regexp(model.genes,'_deleted')))) = false;
        constrainRxn = false(length(rxnInd),1);
        % Figure out if any of the reaction states is changed
        for j = 1:length(rxnInd)
            if (isfield(model, 'rules') && ~isempty(model.rules{rxnInd(j)})) %To avoid errors if the rule is empty
                if (~eval(model.rules{rxnInd(j)}))
                    constrainRxn(j) = true;
                end
            end
        end
        % Constrain flux through the reactions associated with these genes
        if (any(constrainRxn))
            constrRxnNames = model.rxns(rxnInd(constrainRxn));
            if (nargin > 2)
                model = changeRxnBounds(model,constrRxnNames,downRegFraction*model.lb(findRxnIDs(model,constrRxnNames)),'l');
                model = changeRxnBounds(model,constrRxnNames,downRegFraction*model.ub(findRxnIDs(model,constrRxnNames)),'u');
            else
                % Full deletion
                model = changeRxnBounds(model,constrRxnNames,0,'b');
            end
            hasEffect = true;
        end
    end
else
    error(['Gene not in model!']);
end
end