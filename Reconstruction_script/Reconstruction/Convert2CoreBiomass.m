function [model_original] = Convert2CoreBiomass(model_original,strains,inputpath,outputpath)
% This function is to convert models to use the corebiomass based on the
% clade, there should be an input for all phenotype information
% strains should be cell array of strain IDs

fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
    Strain_information(:,i) = data{i};
end
fclose(fid2);
fid2 = fopen('../data/physiology/biomass_type.txt');
format = '%s%s%s%s%s%s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(temp)
    biomass_type(:,i) = temp{i};
end
fclose(fid2);

cd otherchanges/
currentpath = pwd;

% mets to delete from biomass
unusedMets = {'FAD [cytoplasm]','NAD [cytoplasm]','NADH [cytoplasm]','NADP(+) [cytoplasm]',...
    'NADPH [cytoplasm]','riboflavin [cytoplasm]','TDP [cytoplasm]','THF [cytoplasm]','heme a [cytoplasm]','coenzyme A [cytoplasm]'};

% delete those unused mets from panmodel
[~,mets_del_idx] = ismember(unusedMets,model_original.metNames);
bio_rxn_index = find(contains(model_original.rxnNames,'pseudoreaction'));
model_original.S(mets_del_idx,bio_rxn_index) = 0;
fprintf('Panmodel : changed to corebiomass');

for i = 1:length(strains)
    fprintf([strains{i},' : No.',num2str(i),'\n']);
    m = strains{i};
    cd(inputpath)
    reducedModel = load([m,'.mat']);
    cd(currentpath)
    reducedModel = reducedModel.reducedModel;
    
    % this cense is not consense with other field
    if isfield(reducedModel,'csense')
        reducedModel = rmfield(reducedModel,'csense');
    end
    model = reducedModel;
    [~,mets_del_idx] = ismember(unusedMets,model.metNames);
    bio_rxn_index = find(contains(model.rxnNames,'pseudoreaction'));
    model.S(mets_del_idx,bio_rxn_index) = 0;
    
    % rescale biomass according to the proportion
    [~,ID] = ismember(strains{i},Strain_information(:,1));
    type = split(Strain_information(ID,3),',');
    type = type(1);
    [~,ID] = ismember(type,biomass_type(1,:));
    for j = 2:length(biomass_type(:,1))-1 % not input carbonhydrate due to the biomass need to be 1, we will use carbonhydrate to blance out
        model = scaleBioMass(model,biomass_type{j,1},str2double(biomass_type{j,ID}));
    end
    
    % scale out carbohydrate to let biomass to be 1 gram
    [X,~,C,~,~,~,~,~] = sumBioMass(model);
    model = scaleBioMass(model,biomass_type{8,1},(1 - X + C));
    [X,~,C,~,~,~,~,~] = sumBioMass(model);
    cd(outputpath)
    reducedModel = model;
    save([strains{i},'.mat'],'reducedModel')
end
cd(currentpath)
cd ../

end