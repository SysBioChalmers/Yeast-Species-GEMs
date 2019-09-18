function NEWGPRList = AddOrthologGPRrules(model,genes,orthlogs)
%This function converts a gene rule to a logical statement, and then
%update the GPRrules based on orthlog information.
%This function is for reconstrucing the panmodel, called in the function
%UpdatePanGPR.
%       model    A model with GPR rules
%       genes    Genelist which contain genes that has prthlogs
%       orthlogs orthlog information mapping to the genelist
%       rxns     rxns which involived this mapping
%       newGPR   newGPR with adding orthlogs in the GPRrules
%
%
% example:   rxnlist    GPRrules     orthlogs    newGPR
%            rxn1       a or b       a/c         a or b or c
%            rxn2       a and b      a/c         (a and b) or (b and c)
%usage: [rxns,newGPR] = AddOrthologGPRrules(model,genes,orthlogs)
%
%
% Feiran Li 2018.09.25

NEWGPRList = [];
model_r = ravenCobraWrapper(model);
[~,geneindex] =ismember(genes,model_r.genes);
for i = 1:length(geneindex)
    if geneindex(i) ~= 0
        rxns = find(model_r.rxnGeneMat(:,geneindex(i)));
        for j = 1:length(rxns)
            if ~isempty(NEWGPRList) && any(cell2mat(NEWGPRList(:,1))==rxns(j))
                index = find(cell2mat(NEWGPRList(:,1))==rxns(j));
                geneRule = char(NEWGPRList{index,3});
                NewGPR = updateGPRrules(geneRule,genes{i},orthlogs{i});
                NEWGPRList(index,3)=cellstr(NewGPR);
            else
                geneRule = model_r.grRules{rxns(j)};
                NewGPR = updateGPRrules(geneRule,genes{i},orthlogs{i});
                NEWGPRList = [NEWGPRList;rxns(j),cellstr(geneRule),cellstr(NewGPR)];
            end
        end
    end
end

function NewGPR = updateGPRrules(geneRule,gene,orthlog)
%This function converts a gene rule to a logical statement, and then
%asseses if the rule is true (i.e. rxn can still carry flux) or not (cannot
%carry flux).
geneRule = [' ', geneRule, ' '];
geneSets = strsplit(geneRule,'or');
%if ~cellfun(@isempty,strfind(geneSets,gene))
    new = strrep(geneSets,[' ' gene ' '],[' ' orthlog ' ']);
    new = strrep(new,['(' gene ' '],['(' orthlog ' ']);
    new = strrep(new,[' ' gene ')'],[' ' orthlog ')']);
%end
newgpr_temp = unique(union(new,geneSets));
    if length(geneSets) == 1 && length(newgpr_temp) > 1 && ~isempty(strfind(geneRule,'and'))
        NewGPR = strjoin(newgpr_temp,') or (');
        NewGPR = ['(', NewGPR, ')'];
    else
        NewGPR = strjoin(newgpr_temp,'or');
    end
NewGPR = strtrim(NewGPR);
end
end


