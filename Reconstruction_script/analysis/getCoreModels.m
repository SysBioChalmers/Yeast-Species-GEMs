function [Coremodel,VarRxns] = getCoreModels

% CoreRxns
%
%this function is to compare model rxns and mets and genes, and stract core
%models from several models.
%
%   modellist       several models
%
%
%   Coremodel       a structure which contain common rxns and genes and mets
%   Variablerxns    A model structure wich contian variable rxns
%
%   Usage: [Coremodel,Varibalerxns]=CoreRxns(modellist)
%
%   Feiran Li, 2018-09-24

%
cd ..
model = loadYeastModel;
cd ../ComplementaryData/SpecificModelData/
genesMatrix = readtable('genesMatrix_PresenceAbsence_new.xlsx');
StrianData.genes = genesMatrix.geneID;
StrianData.strains = genesMatrix.Properties.VariableNames(2:end)';
StrianData.levels = table2array(genesMatrix(:,2:end));

strain = StrianData.strains;
cd ../../ModelFiles/SSmodels/
rxnexist = zeros(length(model.rxns),length(strain));
geneexist = zeros(length(model.genes),length(strain));
for i = 1:length(strain)
    load([strain{i},'.mat']);
    model1 = reducedModel;
    %[~,index] = ismember(model1.rxns,model.rxns);
    %rxnexist(index(index~=0),i) = 1;
    rxnexist(:,i) = ismember(model.rxns,model1.rxns);
    geneexist(:,i) = ismember(model.genes,model1.genes);
end
CoreRxns = find(all(rxnexist==1,2));
CoreGenes = find(all(geneexist==1,2));
VarRxnsIndex = find(any(rxnexist==0,2));
VarGenesIndex = find(any(geneexist==0,2));
VarRxns = model.rxns(VarRxnsIndex);
VarGenes = model.genes(VarGenesIndex);
Coremodel = removeRxns(model, model.rxns(VarRxnsIndex));
[~,index] = ismember(VarGenes,Coremodel.genes);
Coremodel = removeGenes(Coremodel, Coremodel.genes(index(index~=0)));
cd ../../ComplementaryData/Results
save('VarGenes.mat','VarGenes')
save('VarRxns.mat','VarRxns')
save('rxnexist.mat','rxnexist')
save('geneexist.mat','geneexist')
cd ../../ModelFiles/SSmodels
save('Coremodel.mat','Coremodel')

end




        


