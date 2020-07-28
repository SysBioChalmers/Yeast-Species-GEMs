function acce_orthologs = getAccessGenes(strains,inputpath,type)
% This function is to generate accessory orthologs
% inputpath: the folder of models
% type:      'strict' or 'loose'
%            for when you want to take out genes exist exclusive in accessory
%            reactions

current_path = pwd;
cd(inputpath)
for i = 1:length(strains)
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
[~,Idx] = ismember(acce_rxn,model_original.rxns);
genes_involve = model_original.rxnGeneMat(Idx,:);
genes_out     = model_original.rxnGeneMat(setdiff(1:1:length(model_original.rxns),Idx),:); % core reactions linked with genes 
genes_involve = find(sum(genes_involve,1));% get all access genes index
genes_out = find(sum(genes_out,1));
if strcmp(type,'strict')
    acce_orthologs = model_original.genes(setdiff(genes_involve,genes_out)); % extract the genes only in the access reactions
else
    acce_orthologs = model_original.genes(genes_involve);
end

