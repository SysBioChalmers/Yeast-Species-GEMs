function results = getSubsetTable(indxs,model,solution,varNames,enz)
%get model grRules
if ~enz
    formulas = constructEquations(model,indxs);
    subSystems = {};
    for i=1:length(indxs)
        indx   = indxs(i);
        if isempty(model.subSystems{indx})
            str = ' ';
        else
            str = strjoin(model.subSystems{indx},' // ');
        end
        subSystems = [subSystems; {str}];
    end
    results = table(model.rxns(indxs),model.rxnNames(indxs),formulas,solution(indxs),model.grRules(indxs),subSystems,'VariableNames',varNames);
else
    prots = model.rxnNames(indxs);
    prots = strrep(prots,'prot_','');
    prots = strrep(prots,'draw_','');
    prots = strrep(prots,'_exchange','');
    subSystems = mapEnzymeSubSystems(prots,model);
    results = table(prots,solution(indxs),model.grRules(indxs),subSystems,'VariableNames',varNames);
end
end
     