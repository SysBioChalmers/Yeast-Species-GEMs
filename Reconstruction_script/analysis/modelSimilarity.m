% Figure 2d this function is to get the model similarity analysis for all models we
% generated

inputpath = '../Reconstruction/modelRelated/ssGEMs';
current_path = pwd;
% load group information
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
    Strain_information(:,i) = temp{i};
end
fclose(fid2);
strains = Strain_information(:,1);

% sort as clade
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
group_clade = {};
strains_sortclade = {};
for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    group_clade = [group_clade;Strain_information(idx,2)];
    strains_sortclade = [strains_sortclade;Strain_information(idx,1)];
end

% color code use this one rather the default one in the raven function
color_palette = [166 206 227
106 61 154
31 120 180
178 223 138
51 160 44
251 154 153
227 26 28
253 191 111
255 127 0
202 178 214
254 153 255
177 89 40
0 0 0]/255;

cd(inputpath)
for i = 1:length(strains_sortclade)
    m = strains_sortclade{i};
    load([m,'.mat'])
    reducedModel.id = m;
    models{i} = reducedModel;
end
% compare and generate figure using the function from raven
compStruct = compareMultipleModels(models,true,true,group_clade);
cd(current_path)