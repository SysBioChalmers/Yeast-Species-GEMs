function Essential_cluster

addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/MATLAB-git'))
workdir = pwd;
cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

cd ../../analysis
load('../Reconstruction/StrainData.mat')
strains = StrianData.strains;
nostrain = {'Candida_albicans','Komagataella_pastoris','Saccharomyces_cerevisiae','Schizosaccharomyces_pombe','Yarrowia_lipolytica'}
strains = setdiff(strains,nostrain)
current_path = pwd;
inputpath = '../ModelFiles/mat';

[accurancy,sensitivity,specificity,mcc] = EssentialGeneanalysisnew(strains,inputpath)


end
