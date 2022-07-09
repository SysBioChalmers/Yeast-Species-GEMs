function [pathway,occuranceratio] = findAlternativePWY(strainsquery,target)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to find alternative PWY besides the default pwy in the
% panmodel by serching for all pwys related to the target that you input.
% Step 1: identify all pwys related to the target
% Step 2: load the draft model to check the occurance ratio for each
% pathway
% Step 3: tranform the rxn into the format that the model can add
% Step 4: add those rxns back into the model


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1: identify all pwys related to the target
% identify metacyc rxn database location.
IDfile = 'ravenCobraWrapper.m';
currentPath = pwd;
try
    toolboxPath = which(IDfile);                %full file path
    slashPos    = getSlashPos(toolboxPath);
    toolboxPath = toolboxPath(1:slashPos(end)); %folder path
    cd(toolboxPath);
catch
    disp('raven toolbox cannot be found')
end
cd ../external/metacyc/
load('metaCycRxns.mat') % load all metacyc rxns
cd(currentPath)

% load all metacyc pwy names and IDs
format = '%s %s';
fID       = fopen('../../ComplementaryData/databases/All_pathways_of_MetaCyc.txt');
matrixData  = textscan(fID,format,'Delimiter','\t','HeaderLines',1);
metacyc.pwyIDs      = matrixData{1};
metacyc.pwyNames = matrixData{2};

% find target pwd IDs
pwyidx = contains(metacyc.pwyNames,target);
pwdID = metacyc.pwyIDs(pwyidx);

% Find rxns related to each pathway
total = [];
if ~isempty(pwdID)
    mappedrxns = cell(length(pwdID),1);
    for i = 1:length(pwdID)
        mappedrxns{i} = metaCycRxns.rxns(contains(metaCycRxns.pwys,pwdID(i)));
        total = [total;mappedrxns{i}];
    end
    [a,IA,~] = unique(total);
    sharedrxn = unique(total(setdiff(1:1:length(total),IA)));
end

% Take away shared rxns
if ~isempty(sharedrxn)
    for i = 1:length(pwdID)
        mappedrxns{i} = setdiff(mappedrxns{i},sharedrxn);
    end
end

%% Step 2: load the draft model to check the occurance ratio for each
% pathway
cd ../../draft_GEM_all_yeast_species/
result = [];
cd 'strain specific model from RAVEN_biocyc_55_110'/
for i = 1:length(strainsquery)
    strains = strainsquery(i);
    % load biocyc draft model
    cd(strainsquery{i})
    fid = fopen('draft_GEM.tsv','r');
    draft = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
    draft = [draft{1} draft{2} draft{3} draft{4} draft{5}];
    fclose(fid);
    cd ../
    % index whether the rxn is in the draft
    for j = 1:length(mappedrxns)
        
        [idx,~] = ismember(mappedrxns{j},draft(:,1));
        result = [result;strainsquery{i},pwdID{j},strjoin(mappedrxns{j},';'),num2str(idx),{num2str(sum(idx)/length(idx))}];
    end
end

%% Step 3 




function slashPos = getSlashPos(path)
slashPos = strfind(path,'\');       %Windows
if isempty(slashPos)
    slashPos = strfind(path,'/');   %MAC/Linux
end
end
end
