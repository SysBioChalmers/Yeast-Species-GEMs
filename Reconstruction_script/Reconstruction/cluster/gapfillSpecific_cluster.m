function gapfillSpecific_cluster(a,b)

addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/MATLAB-git'))
setRavenSolver('gurobi')
workdir = pwd;
cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

cd ../
load('aaa.mat')
model_original = panmodel;
load('StrainData.mat')
strains = StrianData.strains(a:b);
inputpath = '../../ModelFiles/mat';
[model_original] = Convert2CoreBiomass(model_original,strains,inputpath,inputpath);
[model_original] = SpecificModelfillGaps(model_original,strains,inputpath);
[model_original] = SubstrateUsageGapFill(model_original,strains,inputpath);
save(['model_origianl',num2str(a),'.mat'],'model_original')
end
