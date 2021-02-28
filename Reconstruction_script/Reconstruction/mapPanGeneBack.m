function [nomatchresult,multimatchresult] = mapPanGeneBack(strains,inputpath,outputpath,changeUniprot)
% This function is to map gene names back to the original genes
if nargin < 4
    changeUniprot = false;
end
currentpath = pwd;
nomatchresult = cell(0,1);
multimatchresult = cell(0,1);
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
    model.genes(idx(idx~=0)) = geneID_core(idx~=0); % fix the geneID by gapfilling so that the strange RNA gene ID can be changed back
    
    para = [];
    for j = 1:length(model.genes)
        a = find(ismember(panID_final,model.genes(j)));
        for m = 1:length(a)
            para = [para;model.genes(j),geneID_core(a(1)),geneID_core(a(m))];
        end
    end   
    
    if ~isempty(para)
    [~,idx] = ismember(para(:,1),model.genes);
    model.genes(idx(idx~=0)) = para((idx~=0),2); % replace gene ID by the first paralog
    para        = para(:,2:3);
    idx         = cellfun(@strcmp,para(:,1),para(:,2));% find paralog 
    para        = para(~idx,:); % generate paralog info for the next step update
    end
    
    cd otherchanges/
    model.grRules = rulesTogrrules(model);
    [grRules,rxnGeneMat]   = standardizeGrRules(model,true);
    model.grRules = grRules;
    model.rxnGeneMat = rxnGeneMat;
    
    % combine the same genes, for example: Sa matched to PANa and PANb
    model = deleteDuplicateGenes(model);
    
    % update the paralog information
    if ~isempty(para)
    model = UpdatePanGPR(para,model); % not update the protein so that panID can be saved there
    end
    
    % find not matched genes and delete those genes, but keep those rxns
    nomatch     = setdiff(model.genes,geneID_core);
    model       = removeGenesFromModel(model,nomatch); % use cobra one so taht the reaction will be kept
 
    %Incorporate a rxnGeneMat consistent with standardized grRules
    model.grRules        = rulesTogrrules(model);
    model = slimGPR(model);
    nomatchresult{i}     = nomatch;
    cd ..
    % get panID as proteins
    model.proteins = model.genes;
    [~,idx] = ismember(geneID_core,model.genes);
    model.proteins(idx(idx~=0)) = panID_final(idx~=0);
    
    
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
    cd(outputpath)
    save([strains{i},'.mat'],'reducedModel')
end

cd(currentpath)
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