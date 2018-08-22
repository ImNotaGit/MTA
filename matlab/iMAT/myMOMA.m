function [solutionDel, totalFluxDiff, constrRxnNames] = myMOMA(WTmodel, v_ref, KOgene)
    % Performs a quadratic version of the MOMA (minimization of metabolic adjustment) approach
    %
    % INPUTS:
    %    WTmodel:            Wild-type model
    %    v_ref:            Reference (i.e wild-type) fluxes
    %    KOgene:           Gene to delete
    %
    % OUTPUTS:
    %    solutionDel:      Deletion solution structure
    %    totalFluxDiff:    Value of the linear MOMA objective, i.e.
    %                      :math:`\sum (v_{wt}-v_{del})^2`
    %    constrRxnNames:   The affected reactions due to gene knock-out
    %
    % myMOMA solves the following:
    %
    % .. math::
    %      min ~&~ \sum (v_{wt} - v_{del})^2 \\
    %          ~&~ S_{del}v_{del} = 0 \\
    %          ~&~ lb_{del} \leq v_{del} \leq ub_{del}
    
    solutionDel.f = [];
    solutionDel.x = [];
    solutionDel.stat = -1;
    totalFluxDiff = [];
    
    % model with KOgene knocked-out
    [model, hasEffect, constrRxnNames] = mydeleteModelGenes(WTmodel, KOgene);
    
    if ~hasEffect
        fprintf('Knocking out the gene does not have any effect on the reactions.\n');
        return
    end
    
    [nMets,nRxns] = size(model.S);
    
    b = zeros(nMets,1);
    A = model.S;
    c = -2*v_ref;
    F = 2*eye(nRxns);
    lb = model.lb;
    ub = model.ub;
    csense(1:nMets) = 'E';
    
    % Solve the linearMOMA problem
    [QPproblem.A,QPproblem.b,QPproblem.F,QPproblem.c,QPproblem.lb,QPproblem.ub,QPproblem.csense,QPproblem.osense] = deal(A,b,F,c,lb,ub,csense,1);
    QPsolution = solveCobraQP(QPproblem, 'printLevel', 0, 'method', 0);
    
    % Get the solution(s)
    if QPsolution.stat == 1
        solutionDel.x = QPsolution.full;
        solutionDel.f = sum(model.c.*solutionDel.x);
        totalFluxDiff = sum((v_ref-solutionDel.x).^2);
    end
    solutionDel.stat = QPsolution.stat;
    solutionDel.solver = QPsolution.solver;
    solutionDel.time = QPsolution.time;

end


function [model, hasEffect, constrRxnNames] = mydeleteModelGenes(model, geneList, downRegFraction)
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
[isInModel,geneNameInd] = ismember(geneList, model.genes_unique_names);

if (isInModel)
  
  geneInd = find(model.genes_unique_map==geneNameInd)
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