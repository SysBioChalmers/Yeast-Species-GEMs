function newModel = getOxPhosGPRs(orgModel)
oxPhos = {'r_1021' 'r_0439' 'r_0438' 'r_0226'};
load('../../models/yeastGEM/Saccharomyces_cerevisiae.mat')
model = ravenCobraWrapper(model);
%Find oxphos rxns in both models
for rxn=oxPhos
    org_Index = find(strcmpi(orgModel.rxns,rxn));
    sce_Index = find(strcmpi(model.rxns,rxn));
    if ~isempty(org_Index)
        grRule = orgModel.grRules{org_Index};
        if isempty(grRule)
            orgModel.grRules{org_Index} = model.grRules{sce_Index};
            genes = model.grRules{sce_Index}
            genes = strrep(genes,' and ',' ');
            genes = strrep(genes,' or ',' ');
            genes = strrep(genes,')','');
            genes = strrep(genes,'(','');
            genes = strsplit(genes,' ')
            length(genes)
            for i = 1:length(genes)
                gene = genes(i);
                if ~contains(orgModel.genes,gene)
                    gene
                    genesToAdd.genes = gene;
                    index = find(strcmpi(model.genes,gene));
                    shortGeneN = model.geneShortNames(index);
                    genesToAdd.geneShortNames = shortGeneN;
                    orgModel = addGenesRaven(orgModel,genesToAdd);
                end
            end
                    
        end
    end
end
newModel = orgModel;
end