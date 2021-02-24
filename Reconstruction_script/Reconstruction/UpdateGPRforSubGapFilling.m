UpdateGPRrules
cd otherchanges/
[a,b,result] = xlsread('../data/geneMining/sub_result.xlsx');
% load old mapping data  geneID with panID
 load('../strainData.mat')
 strains = StrianData.strains;
 % result: mrnaID KO speciesgeneID panID species rxn sub oldGPR

for i = 1:length(strains)    
    fileName    = ['../../../../Multi_scale_evolution/pan_genome/result/id_mapping/',strains{i},'.tsv'];
    fID         = fopen(fileName);
    protData    = textscan(fID,'%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    fclose(fID);
    geneID_core = protData{2}; % geneID in 343 yeast species with @seq
    panID_final = protData{5}; % panID
    draftgeneID      = strcat(protData{6}); % mRNAID geneID CDS
    draftgeneID      = strrep(draftgeneID,' ','&'); % mRNAID&geneID&CDS

    panID_final = strrep(panID_final,'Saccharomyces_cerevisiae@','');
    
    % get the index for genes of the species in the result
    geneIdx = find(endsWith(lower(result(:,2)),lower(strains(i))));
        
    if strcmp(strains{i},'Saccharomyces_cerevisiae') & contains(result(geneIdx,1),'&')
        result(geneIdx,1) = extractBefore(result(geneIdx,1),'&');
    end
    [~,draftGeneIdx] = ismember(result(geneIdx,1), draftgeneID);

    if ~all(draftGeneIdx)
        i
        draftgeneID = extractBefore(draftgeneID,'&');
        [~,draftGeneIdx] = ismember(result(geneIdx,1), draftgeneID);
    end
    result(geneIdx,3) = geneID_core(draftGeneIdx); % fix the geneID by gapfilling so that the strange RNA gene ID can be changed back
    result(geneIdx,4) = panID_final(draftGeneIdx);
    result(geneIdx,5) = repmat(strains(i),length(geneIdx),1); % species
end

% get rxn for ko
fileName    = '/Users/feiranl/Documents/Python/blast/NewGeneMining/rxnlistnew.txt';
fID         = fopen(fileName);
rxn_ko    = textscan(fID,'%s%s%s','Delimiter',',','HeaderLines',1);
fclose(fID);
sub = rxn_ko{1};
rxn = rxn_ko{2};
ko = rxn_ko{3};
ko = strrep(ko,'.fa','');
for i = 1:length(ko)
    ko_tmp = split(ko(i),';');
    for j = 1:length(ko_tmp)
        idx = find(startsWith(result(:,2),ko_tmp));
        result(idx,6) =  rxn(i);
        result(idx,7) =  sub(i);
    end
end
clearvars -except result ko sub rxn model_original
% map the grRules for check and change grRule
strainlst = unique(result(:,5));
load('../../ModelFiles/model_original_new.mat')
inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';
current_path = pwd;
result(:,8) = {''}; %oldgpr
for i = 1:length(strainlst)
    i
    cd(inputpath)
    load([strainlst{i},'.mat'])
    idx = find(ismember(result(:,5),strainlst(i)));
    for j = 1:length(idx)
        rxn_tmp = split(result(idx(j),6),';');
        [~,rxnIdx] = ismember(rxn_tmp,reducedModel.rxns);
        if any(rxnIdx)
            result(idx(j),8) = join(reducedModel.grRules(rxnIdx(rxnIdx~=0)),';');% original grRule
        end
        if isempty(cell2mat(result(idx(j),8)))
            cd(current_path)
            reducedModel = addrxnBack(reducedModel,model_original,rxn_tmp,repmat(result(idx(j),3),length(rxn_tmp),1));
            [~,tmp] = ismember(result(idx(j),3),reducedModel.genes);
            reducedModel.proteins(tmp) = result(idx(j),4);
        end
    end
    cd(inputpath)
    save([strainlst{i},'.mat'],'reducedModel')
end

% remove the rxn without grRule for the rxn listed above
cd(current_path)
rxnlst = unique(result(:,6));
[~,rxnMatrix,~,~,~] = getprecursorMatrixCobra(model_original,strainlst,inputpath,[],0);
removelst = {};
for i = 1:length(rxnlst)
    rxn_tmp = split(rxnlst(i),';');
    idx = find(ismember(result(:,6),rxnlst(i)));
    strain_exclude = setdiff(strainlst,unique(result(idx,5))); % species
    % strains with this rxn
    for j = 1:length(rxn_tmp)
    [~,idx2] = ismember(rxn_tmp(j),model_original.rxns);
    strains_have = strainlst(logical(rxnMatrix(:,idx2)));
    strain_remove = intersect(strains_have,strain_exclude);
    removelst = [removelst;strain_remove,repmat(rxn_tmp(j),length(strain_remove),1)];
    end
end
[~,idx2] = ismember(removelst(:,2),rxn); % index the remove rxn in the ko mapping file
removelst(idx2~=0,3) = sub(idx2(idx2~=0));

tmp = find(idx2 == 0);
idx2 = contains(rxn,removelst(tmp,2));
removelst(tmp,3) = sub(idx2);

[~,~,result_table] = xlsread('../../data/substrateUsageGene.xlsx','RESULTTABLE');
strain_sub_lst = result_table(10:338,1);
result_table = result_table(10:end,2:end);
[~,~,index] = xlsread('../../data/substrateUsageGene.xlsx','index');
[~,idx] = ismember(strainlst,strain_sub_lst); % there are several not in the sublist 332-329
result_table(idx~=0,:) = result_table(idx(idx~=0),:);
rxnexistence = result_table(:,startsWith(index(:,1),'r_'));
FBAresult = result_table(:,endsWith(index(:,1),'_model'));
sub_exp = result_table(:,endsWith(index(:,1),'_exp'));
sub = index(endsWith(index(:,1),'_exp'),1);
sub = extractBefore(sub,'_exp');
sub_exp(strcmp(sub_exp,'n')) = {nan};
sub_exp(strcmp(sub_exp,'v')) = {1};

for i = 1:length(removelst(:,1))
    [~,idx] = ismember(removelst(i,1),strainlst);
    [~,idx2] = ismember(removelst(i,3),sub);
    removelst(i,4) = sub_exp(idx,idx2); % exp result
    removelst(i,5) = FBAresult(idx,idx2); % model prediction
end
change = {};
cd(inputpath)
for i = 1:length(strainlst)
    i
    load([strainlst{i},'.mat'])
    idx = find(ismember(removelst(:,1),strainlst{i}));
    [~,idx2] = ismember(removelst(idx,2),reducedModel.rxns);
    removelst(idx,6) = reducedModel.grRules(idx2);
    model = removeRxns(reducedModel,removelst(idx,2));
    removed = setdiff(reducedModel.rxns,model.rxns); % reducedmodel is old model is new
    change = [change;strainlst(i),join(removed,';'),{'removed'}];
    save([strainlst{i},'.mat'],'reducedModel')
end

% update panmodel
for i = 1:length(rxnlst)
rxn_tmp = split(rxnlst(i),';');
idx = find(ismember(result(:,6),rxnlst(i)));
panGPR = join(unique(result(idx,4)),' or ');
% strains with this rxn
for j = 1:length(rxn_tmp)
[~,idx2] = ismember(rxn_tmp(j),model_original.rxns);
tmp(j) = model_original.grRules(idx2);
if isempty(tmp{j})
    newGPR = panGPR;
else
    newGPR = join([tmp(j),panGPR],' or ');
end
newGPR = sortGPR(newGPR); % fix replicate genes
model_original = addrxnBack(model_original,model_original,rxn_tmp(j),newGPR);
end
y(i,2) = join(tmp,';');
y(i,3) = panGPR;
y(i,4) = join(y(i,2:3),' or ');
clear x tmp
end

grRules = standardizeGrRules(model_original);
model_original.grRules = grRules;
model_original = slimGPR(model_original);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = sortGPR(b)
tempGPR = split(b,' or ');
tempGPR = unique(tempGPR);
tempGPR = setdiff(tempGPR,{''});
a = join(tempGPR,' or ');
end