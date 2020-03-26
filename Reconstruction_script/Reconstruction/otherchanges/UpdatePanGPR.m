function model = UpdatePanGPR(ortholog,model)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UpdatePanGPR
% Add changes from the Pangenome, add orthlogs foe existing genes to update
% GPR in the model
% Input: model, PanGenes.tsv,SGDgeneNames.tsv.
% As for the reference of new GPR, please find detailed information in:
% ComplementaryData/databases/Pangenes.tsv
% NOTE: changeGeneAssociation.m is a function from cobra,
% addOtherlogGenes.m is also called during the function
%
% Feiran Li 2018.09.25
% Feiran Li 2019.09.11 - change to a function ortholog should be supplied
%                         as input format should be a two column cell array : 
%                         ortholog_s288c	panID
%                              YHR092C	    PANXXXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
%Clone the necessary repos:
git('clone https://github.com/SysBioChalmers/yeast-GEM.git')

%Load yeast model:
cd yeast-GEM
model    = load('ModelFiles/mat/yeastGEM.mat');
model    = model.model;
yeastVer = model.modelID(strfind(model.modelID,'_v')+1:end);
cd ..

end
% %change model twice to avoid the useless brackets in the rules
% modelr = ravenCobraWrapper(model);
% model1 = ravenCobraWrapper(modelr);
% model.rules = model1.rules;
% model.grRules = model1.grRules;
%replace the orthologs with the genes that existed in the model to generate
%new GPRs, isoenzymes in the Pan genes for the existed genes in the model
NEWGPRList = AddOrthologGPRrules(model,ortholog(:,1),ortholog(:,2));


for i = 1:length(NEWGPRList(:,1))
    model    = changeGeneAssociation(model, model.rxns{NEWGPRList{i,1}}, NEWGPRList{i,3});
end
% Delete unused genes (if any)
model = removeUnusedGenes(model);


% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = strcat('COBRAProtein',num2str(i));
end

% add gene standard name for new genes
fid = fopen('../../data/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

geneIndex = zeros(1,1);
for i = 1: length(model.genes)
    geneIndex = strcmp(yeast_gene_annotation{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(yeast_gene_annotation{2}{geneIndex})
        model.geneNames{i,1} = yeast_gene_annotation{2}{geneIndex};
    else
        model.geneNames{i,1} = model.genes{i};
    end
end

if isfield(model,'geneShortNames')
    model.geneShortNames = model.genes;
end

%Remove the cloned repos:
%rmdir('yeast-GEM', 's')

% Save model:
%model = rmfield(model,'grRules');
end


