% updateGPRforComplexI
cd otherchanges/
current_path = pwd;
fid2 = fopen('../../data/geneMining/complexI_result.tsv'); % data is generated from geneMining
format = '%s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
result(:,1) = data{1};
result(:,2) = data{1};

% load old mapping data  geneID with panID
 load('../modelRelated/strainData.mat')
 strains = StrianData.strains;
 % result: mrnaID KO_species speciesgeneID panID species rxn sub oldGPR at the end

for i = 1:length(strains)    
    fileName    = ['../Multi_scale_evolution/pan_genome/result/id_mapping/',strains{i},'.tsv']; % this repo is cloned in the building step Multi_scale_evolution, if not exist, then go to clone again to the folder Reconstruction
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

for i = 1:length(strains)
idx = find(ismember(result(:,5),strains(i)));
count(i) = length(idx);
end
strains_withcomplexI = strains(count>=5);

% update panmodel
inputpath = '../modelRelated/ssGEMs';
for i = 1:length(strains_withcomplexI)
cd(inputpath)
load([strains_withcomplexI{i},'.mat']);
idx = find(ismember(result(:,5),strains_withcomplexI(i)));
newGPR = join(result(idx,3),' and ');
cd(current_path)
newGPR = sortGPR(newGPR); % fix replicate genes
reducedModel = addrxnBack(reducedModel,model_original,{'r_5195'},newGPR);
[~,idx2] = ismember(result(idx,3),reducedModel.genes);
reducedModel.proteins(idx2) = result(idx,4);
grRules = standardizeGrRules(reducedModel);
reducedModel.grRules = grRules;
cd(inputpath)
save([strains_withcomplexI{i},'.mat'],'reducedModel')
end
cd(current_path)
cd ../

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = sortGPR(b)
tempGPR = split(b,' or ');
tempGPR = unique(tempGPR);
tempGPR = setdiff(tempGPR,{''});
a = join(tempGPR,' or ');
end