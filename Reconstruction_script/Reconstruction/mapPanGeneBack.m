function [model,nonatchresult] = mapPanGeneBack(strains,inputpath,changeUniprot)
% This function is to map gene names back to the original genes
if nargin < 3
    changeUniprot = false;
end
currentpath = pwd;
nonatchresult = cell(0,1);
for i = 1:length(strains)
    i
    cd(inputpath)
    model = load([strains{i},'.mat']);
    model = model.reducedModel;
    cd(currentpath);
    % load old mapping data  geneID with panID
    fileName    = ['../../../Multi_scale_evolution/pan_genome/result/id_mapping/',strains{i},'.tsv'];
    fID         = fopen(fileName);
    protData    = textscan(fID,'%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    fclose(fID);
    geneID_core = protData{2}; % geneID in 343 yeast species with @seq
    panID_final = protData{5}; % panID
    draftgeneID      = strcat(protData{6},protData{7}); % mRNAID geneID CDS
    draftgeneID      = strrep(draftgeneID,' ','&'); % mRNAID&geneID&CDS
    panID_final = strrep(panID_final,'Saccharomyces_cerevisiae@','');
    
    [~,idx] = ismember(draftgeneID,model.genes);
    model.genes(idx(idx~=0)) = geneID_core(idx~=0); % fix the geneID by gapfilling
    
    para = [];
    for j = 1:length(model.genes)
        a = find(ismember(panID_final,model.genes(j)));
        for m = 1:length(a)
            para = [para;model.genes(j),geneID_core(a(1)),geneID_core(a(m))];
        end
    end
    % replace the model.genes with the first match
    nomatch     = setdiff(model.genes,panID_final);
    model       = removeGenesFromModel(model,nomatch); % use cobra one so taht the reaction will be kept
    
    
    [~,idx] = ismember(para(:,1),model.genes);
    model.genes(idx(idx~=0)) = para((idx~=0),2); % replace gene ID by the first paralog
    para        = para(:,2:3);
    idx         = cellfun(@strcmp,para(:,1),para(:,2));% find paralog 
    para        = para(~idx,:);

    cd otherchanges/
    model.grRules = rulesTogrrules(model);
    [grRules,~]   = standardizeGrRules(model,true);
    model.grRules = grRules;
    
    % update the paralog information
    model = UpdatePanGPR(para,model);
    
    %Incorporate a rxnGeneMat consistent with standardized grRules
    model.grRules        = rulesTogrrules(model);
    [grRules,rxnGeneMat] = standardizeGrRules(model,true);
    model.grRules        = grRules;
    model.rxnGeneMat     = rxnGeneMat;
    nomatchresult{i}     = nomatch;
    cd ..
    if changeUniprot
        fileName       = ['../data/uniprot_geneMap/',strains{i},'.csv'];
        fID            = fopen(fileName);
        protData       = textscan(fID,'%s%s%s','Delimiter','\t','HeaderLines',1);
        SSIDList_core  = protData{1};% geneID in 343 yeast species with @seq
        geneList_final = protData{2};
        fclose(fID);
        
        geneList_final = split(geneList_final,'|');% geneID in uniprot
        geneList_final = geneList_final(:,end); % the last colomn coresponding to be uniprot
 
        [~,idx] = ismember(SSIDList_core,model.genes);
        model.genes(idx(idx~=0)) = geneList_final((idx~=0)); % replace gene ID by the first paralog

        %Incorporate a rxnGeneMat consistent with standardized grRules
        cd otherchanges/
        model.grRules        = rulesTogrrules(model);
        [grRules,rxnGeneMat] = standardizeGrRules(model,true);
        model.grRules        = grRules;
        model.rxnGeneMat     = rxnGeneMat;
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