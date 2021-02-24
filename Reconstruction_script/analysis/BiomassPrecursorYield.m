% This function is to generate the file for decision tree and precursor
% yiled analysis figure 

fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
Strain_information(:,i) = data{i};
end
fclose(fid2);
strains = Strain_information(:,1);

%load panmodel
load('/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/Reconstruction/otherchanges/model_original_withfullgpr.mat')

path = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';
cd ../Reconstruction/otherchanges/
[proMarix,rxnMatrix,mets_test,~,solresult] = getprecursorMatrixCobra(model_original,strains,path);
