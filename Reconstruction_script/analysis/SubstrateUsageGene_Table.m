% this function is to generate a table which compares substrate usage
% result and extract the rxnexistence information
% output is resut_table index
function [result_table,index] = SubstrateUsageGene_Table
current_path = pwd;
inputpath = '../Reconstruction/modelRelated/ssGEMs';

% read the index for the resulttable
[~,~,index] = xlsread('../Reconstruction/modelRelated/substrateUsageGene.xlsx','index');

% read the model and data
load('../Reconstruction/modelRelated/panModel.mat');

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
SubModelName(1:7) = strcat(SubModelName(1:7),'_fermentation');
SubCondition = data(:,3);
Subtype = data(:,4);
data(:,1:4) = [];
fclose(fid2);

% generate the result
[~,FBAresult] = SubstrateUsageTest(model_original,false,inputpath);
cd ../Reconstruction/otherchanges/
[~,rxnMatrix,~,~,~] = getprecursorMatrixCobra(model_original,strainlist,inputpath,[],0);


% generate the table as index said
result_table = cell(length(strainlist),length(index(:,1)));
KEGG = cell(1,length(index(:,1)));
MNX = cell(1,length(index(:,1)));
for i = 1:length(index(:,1))
    i
    if endsWith(index{i,1},'_model')
        trait = replace(index{i,1},'_model','');
        [~,idx] = ismember(trait,SubModelName);
        result_table(:,i) = num2cell(FBAresult(idx,:)');
    elseif endsWith(index{i,1},'_exp')
        trait = replace(index{i,1},'_exp','');
        [~,idx] = ismember(trait,SubModelName);
%         temp = data(idx,:)';
%         temp = strrep(temp,'v','1');
%         temp = strrep(temp,'n','nan');
%         temp = cellfun(@str2num, temp, 'UniformOutput', false);
        result_table(:,i) = data(idx,:)';
    else
        rxn = split(index{i,1},';');
        [~,idx] = ismember(rxn,model_original.rxns);
        if idx~=0
            existence = rxnMatrix(:,idx);
            result_table(:,i) = num2cell(sum(existence,2));
            KEGG(1,i) = join(model_original.rxnKEGGID(idx),';');
            MNX(1,i) = join(model_original.rxnMetaNetXID(idx),';');
        end
    end
end

% update the model info  remind should take care of all zeros for rxn
pfisher = nan([1,length(index(:,1))]);
for i = 1:length(index(:,1))
if endsWith(index(i,1),'_exp')
X_exp = result_table(:,i);
X_exp = strrep(X_exp,'n','nan');
X_exp = strrep(X_exp,'v','1');
X_exp = cell2mat(cellfun(@str2num, X_exp, 'UniformOutput', false));
elseif  startsWith(index(i,1),'r_') || endsWith(index(i,1),'_model')
X_model = cell2mat(result_table(:,i));
X_model(abs(X_model) > 0) = 1;
% X2 claculation and get p from fisher test
tp(i) = length(X_exp(X_exp>0 & X_model>0));
tn(i) = length(X_exp(X_exp==0 & X_model==0));
fp(i) = length(X_exp(X_exp==0 & X_model>0));
fn(i) = length(X_exp(X_exp>0 & X_model==0));
[tbl,x2(i),p(i),lables] = crosstab(X_exp,X_model);
if numel(tbl) == 4
    [h(i),pfisher(i),~] = fishertest(tbl);
elseif numel(tbl) == 2
    tbl = [tbl,[0;0]]; % model rxn all 1
    [h(i),pfisher(i),~] = fishertest(tbl);
else
    tbl = [[tbl;0],[0;0]]; % glc
    [h(i),pfisher(i),~] = fishertest(tbl);
end
end
end

accuracy = (tp + tn)./(fn + fp + tp + tn);
sensitivity = (tp./(tp+fn));
specificity = (tn./(tn+fp));
positivePredictive = (tp./(tp+fp));
negativePredictive = (tn./(fn+tn));
index = [index(:,1),num2cell([pfisher',accuracy',sensitivity'])];

cd(current_path)