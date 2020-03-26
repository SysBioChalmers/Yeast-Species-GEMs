function [model,noexist,ratio] = mapPanGeneBack(strains,changeUniprot,inputpath)
% This function is to map gene names back to the original genes
if nargin < 2
    changeUniprot = false;
end
currentpath = pwd;

for i = 1:length(strains) 
    i
    % load old mapping data  geneID with panID
    fileName = ['../../../Multi_scale_evolution/pan_genome/result/id_mapping/',strains{i},'.tsv'];
    fID       = fopen(fileName);
    protData  = textscan(fID,'%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    geneID_core      = protData{2}; % geneID in 343 yeast species with @seq
    panID_final     = protData{5}; % panID 
    panID_final = strrep(panID_final,'Saccharomyces_cerevisiae@','');
    
    m = strains{i};
    cd(inputpath);
    model = load([m,'.mat']);
    cd(currentpath);
    model = model.reducedModel;
    para = [];
    for j = 1:length(model.genes)
        Idx = find(strcmp(model.genes(j),panID_final));
        if length(Idx) > 1
            model.genes(j) = geneID_core(Idx(1));
            para = [para;repmat(geneID_core(Idx(1)),length(Idx)-1,1),geneID_core(Idx(2:end),:)];
        else
            model.genes(j) = geneID_core(Idx(1));
        end 
    end
    
    [~,Idx] = ismember(model.genes,panID_final);
    noexist = setdiff(model.genes,panID_final);
    model.genes(Idx~=0) = geneID_core(Idx(Idx~=0));
    model.grRules=rulesTogrrules(model);
    [grRules,rxnGeneMat] = standardizeGrRules(model,true);
    model.grRules      = grRules;
    
    % update the paralog information
    model = UpdatePanGPR(para,model);
     %Incorporate a rxnGeneMat consistent with standardized grRules
     model.rxnGeneMat   = rxnGeneMat;
     ratio = length(Idx(Idx~=0))/length(model.genes);
     
    if changeUniprot
        fileName = ['../data/geneMap/',strains{i},'.csv'];
        fID       = fopen(fileName);
        protData  = textscan(fID,'%s%s%s','Delimiter','\t','HeaderLines',1);
        SSIDList_core      = protData{1};% geneID in 343 yeast species with @seq
        geneList_final     = protData{2};
        geneList_final = split(geneList_final,'|');% geneID in uniprot
        geneList_final = geneList_final(:,3); % the last colomn coresponding to be uniprot
        [~,Idx] = ismember(model.genes,SSIDList_core);
        noexist = setdiff(model.genes,SSIDList_core);
        model.genes(Idx~=0) = geneList_final(Idx(Idx~=0));
        model.grRules=rulesTogrrules(model);
        [grRules,rxnGeneMat] = standardizeGrRules(model,true);
        model.grRules      = grRules;
        %Incorporate a rxnGeneMat consistent with standardized grRules
        model.rxnGeneMat   = rxnGeneMat;
        ratio = length(Idx(Idx~=0))/length(model.genes);
    end   
    reducedModel = model;
    save([strains{i},'.mat'],'reducedModel')
end


function grRules=rulesTogrrules(model)
%This function takes rules, replaces &/| for and/or, replaces the x(i)
%format with the actual gene ID, and takes out extra whitespace and
%redundant parenthesis introduced by COBRA, to create grRules.
grRules = strrep(model.rules,'&','and');
grRules = strrep(grRules,'|','or');
for k = 1:length(model.genes)
    grRules = strrep(grRules,['x(' num2str(k) ')'],model.genes{k});
end
grRules = strrep(grRules,'( ','(');
grRules = strrep(grRules,' )',')');
grRules = regexprep(grRules,'^(',''); %rules that start with a "("
grRules = regexprep(grRules,')$',''); %rules that end with a ")"
end
end