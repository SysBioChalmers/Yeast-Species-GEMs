% this function is to generate a table which compares substrate usage
% result and extract the rxnexistence information

inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';

% read the index for the resulttable
[~,~,index] = xlsread('../data/substrateUsageGene1.xlsx','index');

% read the model and data
load('../Reconstruction_script/ModelFiles/model_original_new.mat')

fid2 = fopen('../data/physiology/Biolog_substrate.tsv');
format = repmat('%s ',1,333);
format = strtrim(format);
substrate = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(substrate)
    data(:,i) = substrate{i};
end
strainlist = data(1,5:end);
data(1,:) = [];
SubBiologName = data(:,1);
SubModelName = data(:,2);
SubCondition = data(:,3);
Subtype = data(:,4);
data(:,1:4) = [];
fclose(fid2);

% generate the result
[accuracy,FBAresult] = SubstrateUsageTest(model_original,false,inputpath);
cd ../Reconstruction/otherchanges/
[~,rxnMatrix,~,~,~] = getprecursorMatrixCobra(model_original,strainlist,inputpath,[],0);


% generate the table as index said
result_table = cell(length(strainlist),length(index(:,1)));
for i = 1:length(index(:,1))
    i
    if endsWith(index{i,1},'_model')
        trait = replace(index{i,1},'_model','');
        [~,idx] = ismember(trait,SubModelName);
        result_table(:,i) = num2cell(FBAresult(idx,:)');
    elseif endsWith(index{i},'_exp')
        trait = replace(index{i,1},'_exp','');
        [~,idx] = ismember(trait,SubModelName);
%         temp = data(idx,:)';
%         temp = strrep(temp,'v','1');
%         temp = strrep(temp,'n','nan');
%         temp = cellfun(@str2num, temp, 'UniformOutput', false);
        result_table(:,i) = data(idx,:)';
    else
        rxn = split(index{i},';');
        [~,idx] = ismember(rxn,model_original.rxns);
        if idx~=0
            existence = rxnMatrix(:,idx);
            result_table(:,i) = num2cell(sum(existence,2));
        end
    end
end


