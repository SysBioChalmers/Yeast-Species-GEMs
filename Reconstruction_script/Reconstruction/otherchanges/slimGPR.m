function model = slimGPR(model)
%This function slim the grRules of duplicate GPRs
%This function is for reconstrucing the panmodel, after updating all
%paralog information in gprs, we should check the gpr duplication
%       model    A model with GPR rules
%
%
% example:   rxnlist    GPRrules         newGPR
%            rxn1       a or b or b     a or b
%usage: model = slimGPR(model)
%
%
% Feiran Li 2020.03.01


if ~isfield(model,'grRules') 
    model = ravenCobraWrapper(model);
end

for i = 1:length(model.grRules)
    GPR = model.grRules{i};
    if ~isempty(GPR)
        genesSets = getSimpleGeneSets(GPR);
        GeneSetMatrix = zeros(length(genesSets),length(model.genes));
        for j = 1:length(genesSets)
            subunit  = transpose(strsplit(genesSets{j},' and '));
            subunit = strtrim(subunit);
            for k = 1:length(subunit)
                genePos = find(strcmpi(model.genes,subunit(k)));
                if ~isempty(genePos)
                    GeneSetMatrix(j,genePos) = 1;
                end
            end
        end
        [~,ia,~] = unique(GeneSetMatrix,'rows');
        genesSets = genesSets(sort(ia));
        if length(genesSets) == 1
            NewGPR = cell2mat(genesSets);
        else
            NewGPR = strjoin(genesSets,') or (');
            NewGPR = ['(', NewGPR, ')'];
            NewGPR = strrep(NewGPR,'  ',' ');
        end
        NewGPR = strtrim(NewGPR);
        model.grRules{i} = NewGPR;
    end
   
end

[grRules,rxnGeneMat] = standardizeGrRules(model,false); % change '(A) or (B)' to 'A or B'
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
rules = grrulesToRules(model);
model.rules = rules;
model = buildRxnGeneMat(model);
model = removeUnusedGenes(model);


%Function that gets a cell array with all the simple geneSets in a given
%grRule string
function genesSets = getSimpleGeneSets(originalSTR)
genesSets  = [];
%If gene rule is not empty split in all its different isoenzymes
if ~isempty(originalSTR)
    originalSTR = strtrim(originalSTR);
    %Remove all brackets
    originalSTR = [' ', originalSTR, ' '];
    %Remove all brackets
    originalSTR = strrep(originalSTR,'(','');
    originalSTR = strrep(originalSTR,')','');
    originalSTR = strrep(originalSTR,' or ',' &% ');
    %Split all the different genesSets
    genesSets = strsplit(originalSTR,'&%');
    genesSets = unique(genesSets,'stable');
end
end

function rules=grrulesToRules(model)
%This function just takes grRules, changes all gene names to
%'x(geneNumber)' and also changes 'or' and 'and' relations to corresponding
%symbols
replacingGenes=cell([size(model.genes,1) 1]);
for m=1:numel(replacingGenes)
    replacingGenes{m}=strcat('x(',num2str(m),')');
end
rules = strcat({' '},model.grRules,{' '});
for m=1:length(model.genes)
    rules=regexprep(rules,[' ' model.genes{m} ' '],[' ' replacingGenes{m} ' ']);
    rules=regexprep(rules,['(' model.genes{m} ' '],['(' replacingGenes{m} ' ']);
    rules=regexprep(rules,[' ' model.genes{m} ')'],[' ' replacingGenes{m} ')']);
end
rules=regexprep(rules,' and ',' & ');
rules=regexprep(rules,' or ',' | ');
rules=strtrim(rules);
rules = strrep(rules,'  ',' ');
end
end