function model = deleteDuplicateGenes(model) 
% this function is to delete duplicated genes and combine rules
    [uniquegene, ind] = unique(model.genes, 'stable');
    % duplicate indices
    duplicate_ind = setdiff(1:length(model.genes), ind);
    if ~isempty(duplicate_ind)
        [~,tmp] = ismember(model.genes(duplicate_ind),uniquegene);
        new_ind = ind(tmp);
        for m = 1:length(duplicate_ind)
            rxnidx = find(contains(model.rules,['x(',num2str(duplicate_ind(m)),')']));
            model.rules(rxnidx(rxnidx~=0)) = strrep(model.rules(rxnidx(rxnidx~=0)),['x(',num2str(duplicate_ind(m)),')'],['x(',num2str(new_ind(m)),')']);
            model.genes(duplicate_ind(m)) = {[model.genes{duplicate_ind(m)},'deleted']};
        end

        % fix other field
        model.grRules        = rulesTogrrules(model);
        [grRules,rxnGeneMat] = standardizeGrRules(model,true);
        model.grRules        = grRules;
        model.rxnGeneMat     = rxnGeneMat;
        model = removeUnusedGenes(model);
    else
        warning('there is no duplicate genes exist')
    end