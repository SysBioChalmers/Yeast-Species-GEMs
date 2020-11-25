function getmodel_cluster(a,b)
% a b is the strain index that you want to 
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
workdir = pwd;
cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

cd ../
load('getmodel_structureforcluster.mat')
cd otherchanges/
for i = a:b
    disp([StrianData.strains{i},' No.',num2str(i)])
    [~,ID] = ismember(StrianData.strains(i),StrianData.strains);
    lvl = StrianData.levels(:,ID);
    lvl_temp = logical(lvl);
    genesExist = StrianData.genes(lvl_temp);
    [~,ID] = ismember(genesExist,ortholog(:,2));
    ortholog_strian = ortholog(ID(ID~=0),:);
    model_temp = model;
    if ~isempty(ortholog_strian)
    model_temp = UpdatePanGPR(ortholog_strian,model);
    end
    [reducedModel,resultfile] = SpecificModel(model_temp,StrianData,StrianData.strains(i),'../../ModelFiles/mat');
end