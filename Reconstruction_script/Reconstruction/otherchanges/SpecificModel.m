function [reducedModel,resultfile] = SpecificModel(model,StrianData,strain,inputpath)
% This function is to generate strain specific models from a panmodel.
% Input of this function are: a panmodel, which will be loaded using the
% function as loadYeastModel and a gene_strainMatrix, which will be loaded
% from ../ComplementaryData/genesMatrix_PresenceAbsence_new.xlsx.
%Usage: [reducedModel,resultfile] = SpecificModel(strain)
%       of [reducedModel,resultfile] = SpecificModel  %note this one
%       generate SSmodels for 343 yeast/fungi species

% load presenceAvsence data
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

% If the supplied object is a character array, then convert it to a cell array
if ischar(strain)
    strain={strain};
end

% Check that the strain exists
if ~ismember(upper(strain),upper(StrianData.strains))
    EM='The strain name does not match in the inputlist';
    dispEM(EM);
end

currentpath = pwd;
resultfile = [];

%% generate the genelist that don't exist
for j = 1:length(strain)
    
    [~,ID] = ismember(strain(j),StrianData.strains);
    lvl = StrianData.levels(:,ID);
    lvl_tmp = lvl == 0;
    idx = ismember(upper(StrianData.genes),upper(model.genes));
    genelist = StrianData.genes(lvl_tmp & idx);
     
    %generate the specific model according to type
    [reducedModel,resultfile] = removeGenesWithRatio(model,genelist,true,true,false,0.5);%add annotaion and generate result file
    reducedModel.id=[strain{j},' specific model genereted from panYeast'];
    
    % change the ID to it specific IDs
    if isfield(reducedModel,'metSBOTerms')
    reducedModel = rmfield(reducedModel,'metSBOTerms');
    reducedModel = rmfield(reducedModel,'rxnSBOTerms');
    end
    cd(inputpath)
    save([strain{j},'.mat'],'reducedModel')
    save([strain{j},'result.mat'],'resultfile')
    cd(currentpath)
end
end
