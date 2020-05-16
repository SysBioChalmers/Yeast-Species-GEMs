function model = AddMissingOrthologsInYeastGEM(model,mapping)
%This function is to change genes that only exist in GEM to PanID. There are only
%one geneID kept in PanID for duplicate genes.
% This is to fill the gap for the yeastGEM, actually these oethologs should
% exist in yeastGEM but not. Later if this gap is filled in yeastGEM, we
% can skip this step.
%input: a new gpr list for rxns: Missingortholog_S288c, this is generated based on GenesOnlyinGEMPanIDMapping.tsv; a cobra model
%output: new model with updated GPR rules.

if nargin < 2
%Load mapping list:
fid = fopen('../../ComplementaryData/SpecificModelData/Missingortholog_S288c.tsv');
mapplist = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
mapping.rxnIDs  = mapplist{3};
mapping.new_GPR  = mapplist{5};
fclose(fid);
end

[~,rxnindex] = ismember(mapping.rxnIDs,model.rxns);
for i = 1:length(rxnindex)
    model = changeGeneAssociation(model, model.rxns{rxnindex(i)},mapping.new_GPR{i});
end

% Add gene standard name for new genes
fid = fopen('../../data/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);
for i = 1: length(model.genes)
    geneIndex = strcmp(yeast_gene_annotation{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(yeast_gene_annotation{2}{geneIndex})
        model.geneNames{i} = yeast_gene_annotation{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end

% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = strcat('COBRAProtein',num2str(i));
end


% Save model:
model = removeUnusedGenes(model);
%cd ../
%saveYeastModel(model)

end
