% This function is to generate the file for decision tree and precursor
% yiled analysis 
currentpath = pwd;
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
Strain_information(:,i) = data{i};
end
fclose(fid2);
strains = Strain_information(:,1);

%load panmodel
load('../modelRelated/panModel.mat')

path = '../modelRelated/ssGEMs';
cd ../Reconstruction/otherchanges/
[proMarix,rxnMatrix,mets_test,~,solresult] = getprecursorMatrixCobra(model_original,strains,path);
matrix.proMarix = proMarix;
matrix.rxnMatrix = rxnMatrix;
matrix.species = strains;
matrix.mets_test = mets_test;
save('../modelRelated/allMatrix.mat','matrix')
cd(currentpath)
