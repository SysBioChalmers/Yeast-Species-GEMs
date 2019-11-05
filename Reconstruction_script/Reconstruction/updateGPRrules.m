function NewGPR = updateGPRrules(geneRule,gene,orthlog)
%This function converts a gene rule to a logical statement, and then
%asseses if the rule is true (i.e. rxn can still carry flux) or not (cannot
%carry flux).
geneRule = [' ', geneRule, ' '];
geneSets = strrep(geneRule,[' or '],[' &% ']);
geneSets = strsplit(geneSets,'&%');
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
        newgpr_temp = strtrim(newgpr_temp);
        NewGPR = strjoin(newgpr_temp,' or ');
    end
NewGPR = strtrim(NewGPR);
end