function savemodel_cluster(a,b)

addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/MATLAB-git'))
workdir = pwd;
cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

cd ../
load('StrainData.mat')
strains = StrianData.strains(a:b);
current_path = pwd;
inputpath = '/cephyr/users/feiranl/Hebbe/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat/'
%[nomatchresult,multimatchresult] = mapPanGeneBack_cluster(strains,inputpath,inputpath,false);
cd(current_path);
for i = 1:length(strains)
disp(['strains: ',num2str(i+a-1),strains{i}])
cd(inputpath)
load([strains{i},'.mat'])
cd(current_path);
cd ../
saveSSModel(reducedModel)
end

end
