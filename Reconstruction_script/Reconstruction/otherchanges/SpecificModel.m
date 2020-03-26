function [reducedModel,resultfile] = SpecificModel(model,StrianData,strain,inputpath)
% This function is to generate strain specific models from a panmodel.
% Input of this function are: a panmodel, which will be loaded using the
% function as loadYeastModel and a gene_strainMatrix, which will be loaded
% from ../ComplementaryData/genesMatrix_PresenceAbsence_new.xlsx.
%Usage: [reducedModel,resultfile] = SpecificModel(strain)
%       of [reducedModel,resultfile] = SpecificModel  %note this one
%       generate SSmodels for 1011 strians
% %load model
% cd ..
% model = loadYeastModel;
%change model to raven format
% if isfield(model,'grRules')
% model = rmfield(model,'grRules');
% end
% if isfield(model,'rxnGeneMat')
% model = rmfield(model,'rxnGeneMat');
% end
% modelr = ravenCobraWrapper(model);
% model = modelr;



%load presenceAvsence data
if nargin <2
genesMatrix = readtable('../ComplementaryData/SpecificModelData/genesMatrix_PresenceAbsence_new.xlsx');
StrianData.genes = genesMatrix.geneID;
StrianData.strains = genesMatrix.Properties.VariableNames(2:end)';
StrianData.levels = table2array(genesMatrix(:,2:end));
end

if nargin<3
    strain = StrianData.strains;
end
if nargin < 4
    inputpath = pwd;
end
%If the supplied object is a character array, then convert it to a cell
%array
if ischar(strain)
    strain={strain};
end

%Check that the strain exists
if ~ismember(upper(strain),upper(StrianData.strains))
    EM='The strain name does not match';
    dispEM(EM);
end

currentpath = pwd;
resultfile = [];

% %generate the genelist that don't exist
for j = 1:length(strain)
    
    [~,ID] = ismember(strain(j),StrianData.strains);
    lvl = StrianData.levels(:,ID);
    lvl_tmp = lvl == 0;
    idx = ismember(upper(StrianData.genes),upper(model.genes));
    genelist = StrianData.genes(lvl_tmp & idx);



%generate the specific model according to type
reducedModel = removeGenes(model,genelist,true,true,false);%add annotaion and generate result file
reducedModel.id=[strain{j},' specific model genereted from panYeast'];

% change the ID to it specific IDs



cd(inputpath)
save([strain{j},'.mat'],'reducedModel')
cd(currentpath)
resultfile = [resultfile;strain(j),length(genelist),length(reducedModel.genes),length(reducedModel.rxns),length(reducedModel.mets)];
end

%cd ..
fid2 = fopen('specificModelResultFile.tsv','w');
formatSpec = '%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fid2,formatSpec,'strain','genelistnumber','genes','rxns','mets','sol.x');
for i = 1:length(resultfile(:,1))
    fprintf(fid2,formatSpec,char(resultfile(i,1)),num2str(resultfile{i,2}),num2str(resultfile{i,3}),num2str(resultfile{i,4}),num2str(resultfile{i,5}));
end
fclose(fid2);
end
