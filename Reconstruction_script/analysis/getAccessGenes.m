function genes_ortholog = getAccessGenes(strains,inputpath,type)

current_path = pwd;
cd(inputpath)
for i = 1:343
    m = strains{i};
    load([strains{i},'.mat'])
    model = reducedModel;
    if i == 1
        core_rxn = model.rxns;
        acce_rxn = model.rxns;
        pan_rxn = model.rxns;
    else
        core_rxn = intersect(core_rxn,model.rxns);
        pan_rxn = union(model.rxns,pan_rxn);
        acce_rxn = setdiff(pan_rxn,core_rxn);
    end
end
acce_rxn(2251) = {'r_5225'}; % lysine degradation
acce_rxn(2039) = []; % beta_cello
[~,Idx] = ismember(acce_rxn,model_original.rxns);
genes_involve = model_original.rxnGeneMat(Idx,:);
genes_out     = model_original.rxnGeneMat(setdiff([1:1:length(model_original.rxns)],Idx),:);
genes_involve = find(sum(genes_involve,1));
genes_out = find(sum(genes_out,1));
if strcmp(type,'strict')
    genes_ortholog = model_original.genes(setdiff(genes_involve,genes_out));
else
    genes_ortholog = model_original.genes(genes_involve);
end

