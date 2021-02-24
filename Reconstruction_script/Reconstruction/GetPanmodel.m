% This function is the main function for generating panmodel and specific
% models
%% Preparation
addpath(genpath('/Users/feiranl/Documents/GitHub/RAVEN'))
addpath(genpath('/Users/feiranl/Documents/GitHub/cobratoolbox'))
addpath(genpath('/Users/feiranl/Documents/GitHub/MATLAB-git'))
addpath(genpath('/Users/feiranl/Library/gurobi901/mac64/matlab/'))

initCobraToolbox
git('clone https://github.com/SysBioChalmers/yeast-GEM.git')

%Load yeast model:
cd yeast-GEM
model    = load('ModelFiles/mat/yeastGEM.mat');
model    = model.model;
cd ..

% modify the rules field to avoid the useless brackets
model1 = ravenCobraWrapper(model);
grRules = model1.grRules;
model1 = ravenCobraWrapper(model1);
model.rules = model1.rules;
model.grRules = grRules; %get a field of grRules
clear model1

% expand the grRules to be paralogs in each gpr rules
git('clone https://github.com/SysBioChalmers/Multi_scale_evolution.git')

fileName = 'Multi_scale_evolution/pan_genome/result/id_mapping/Saccharomyces_cerevisiae.tsv';
fID       = fopen(fileName);
protData  = textscan(fID,'%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
geneID_core      = protData{2}; % geneID in 343 yeast species with @seq
panID_final     = protData{5}; % panID
panID_final = strrep(panID_final,'Saccharomyces_cerevisiae@','');
geneID_core = strrep(geneID_core,'Saccharomyces_cerevisiae@','');
[~,ID] = ismember(geneID_core,model.genes);
para = [model.genes(ID(ID~=0)),panID_final(ID~=0)];
para_idx = setdiff(para(:,1),para(:,2)); % find the paralogs in sce
[~,ID] = ismember(para_idx,para(:,1));
para = para(ID,:);

for i = 1:length(para(:,1))
    if ~ismember(para(i,2),model.genes)
        [~,ID] = ismember(para(i,1),model.genes);
        model.genes(ID) = para(i,2);
    end
    model.grRules = strrep(model.grRules,[para{i,1},' '],[para{i,2},' ']);
    model.grRules = strrep(model.grRules,[para{i,1},')'],[para{i,2},')']);
    ID = find(endsWith(model.grRules,para{i,1}));
    model.grRules(ID(ID~=0)) = strrep(model.grRules(ID(ID~=0)),para{i,1},para{i,2});
end
%cd Yeast-Species-GEMs/Reconstruction_script/Reconstruction/
cd otherchanges/
model = slimGPR(model); % now the model.grRules has been changed to only the representative IDs
cd ../

%change orthlog information to the representative IDs
fid      = fopen('../../find_homolog_for_panID_py/result/pan_hit_mapping_panYeast_v2_PI@70.tsv');
orth     = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
ortholog(:,1)     = orth{1};
ortholog(:,2) = orth{2};
for i = 1:length(para(:,1))
    idx = ismember(ortholog(:,1),para(i,1));
    ortholog(idx,1) = para(i,2); %replace the ortholog with panID
end
fclose(fid);
save('ortholog_changedtoPanID.mat','ortholog');

clearvars -except model ortholog geneID_core panID_final
%% Panmodel expansion
% Three input: new rxn; ortholog information for existing rxns; strain x
% gene matrix.
% load new rxn and new metabolites and generate three tsv files for next step: adding new rxns and mets into the model
% mapping metaNetIDs
format = '%s %s %s %s %s %s %s %s %s %s ';
fID       = fopen('../../rxn_annotate/result/new_met_information_from MNX_database.txt');
matrixData  = textscan(fID,format,'Delimiter','\t','HeaderLines',1);
matrix.rxnIDs      = matrixData{1};
matrix.mettype = matrixData{2};
matrix.metcoef  = matrixData{3};
matrix.metcompartments = repmat({'cytoplasm'},length(matrix.metcoef),1);
newmet.metNames         = matrixData{7};
newmet.metFormulas      = matrixData{5};
newmet.metCharges       = cellfun(@str2num,replace(matrixData{6},'NA','0'));
newmet.metKEGGID        = matrixData{9};
newmet.metChEBIID       = matrixData{8};
newmet.metMetaNetXID    = matrixData{4};
fclose(fID);

% Matching newmat with existing mets in the model through mapping MNXID,
% CHEBI ID and KEGG ID.
[~,ID] = ismember(newmet.metMetaNetXID,model.metMetaNetXID);
metname_temp = split(model.metNames(ID(ID~=0)),' [');
metname_temp = cellstr(metname_temp(:,1));
newmet.metNames(ID~=0) = metname_temp;
matrix.metIDs = newmet.metNames;

fID       = fopen('../../rxn_annotate/result/new_rxn_information_from MNX_database.txt');
rxnData = textscan(fID,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newrxn.ID  = rxnData{1};
newrxn.Rev = cellfun(@str2num,replace(rxnData{24},'need_manual_check','1'));
% newrxn.GPR = replace(rxnData{8},';',' or ');
% newrxn.GPR = replace(newrxn.GPR,'NA',',');
% newrxn.GPR = replace(newrxn.GPR,'"','');
newrxn.GPR(:,4) = replace(rxnData{6},'"',''); %'panID_RAVEN_biocyc'
newrxn.GPR(:,3) = replace(rxnData{5},'"',''); %'panID_RAVEN_kegg'
newrxn.GPR(:,2) = replace(rxnData{4},'"',''); %'panID_eggnog_web'
newrxn.GPR(:,1) = replace(rxnData{3},'"',''); %'panID_kegg_web'
for i = 1:length(newrxn.ID)
tmp = join(newrxn.GPR(i,1:3),';');
tempGPR = split(tmp,';');
tempGPR = unique(tempGPR);
tempGPR = setdiff(tempGPR,{''});
tempGPR = setdiff(tempGPR,{'NA'});
tempGPR = join(tempGPR,' or ');
Note(i,1) = 1;
if strcmp(tempGPR,'')|| isempty(tempGPR)
tmp = newrxn.GPR(i,4);
tempGPR = split(tmp,';');
tempGPR = unique(tempGPR);
tempGPR = setdiff(tempGPR,{''});
tempGPR = setdiff(tempGPR,{'NA'});
Note(i,1) = 2;
tempGPR = join(tempGPR,' or ');
end
GPR(i,1) = tempGPR;
end
newrxn.GPR = GPR;
newrxn.GPRNote = Note;
newrxn.rxnNames     = rxnData{12};
newrxn.rxnNamesKEGG     = rxnData{11};
newrxn.rxnpathway = rxnData{27}; % please do make sure that happen
% find empty rxnnames, 1) replace that with keggrxnname 2) replace that
% with MNXID
idx=cellfun('isempty',newrxn.rxnNames);
newrxn.rxnNames(idx) = newrxn.rxnNamesKEGG(idx);
newrxn.rxnNames = replace(newrxn.rxnNames,'NA','');
idx=cellfun('isempty',newrxn.rxnNames);
newrxn.rxnNames(idx) = newrxn.ID(idx);
newrxn.rxnECNumbers = replace(rxnData{16},'NA','');
newrxn.rxnKEGGID    = replace(rxnData{10},'kegg:','');
newrxn.rxnKEGGID    = replace(newrxn.rxnKEGGID,'NA','');
newrxn.rxnpathway   = replace(newrxn.rxnpathway,'NA','');
newrxn.rxnMetaNetXID   = rxnData{1};
fclose(fID);

% fix the NA in coef/charge
NAidx = find(contains(matrix.metcoef,'n'));
NArxn = unique(matrix.rxnIDs(NAidx));
NArxnidx = find(contains(matrix.rxnIDs,NArxn));
matrix.rxnIDs(NArxnidx) = [];
matrix.mettype(NArxnidx) = [];
matrix.metcoef(NArxnidx)  = [];
matrix.metIDs(NArxnidx)  = [];
matrix.metcoef = cellfun(@str2num,matrix.metcoef);
matrix.metcompartments(NArxnidx) = [];
newmet.metNames(NArxnidx)         = [];
newmet.metFormulas(NArxnidx)      = [];
newmet.metCharges(NArxnidx)       = [];
newmet.metKEGGID(NArxnidx)        = [];
newmet.metChEBIID(NArxnidx)       = [];
newmet.metMetaNetXID(NArxnidx)    = [];

NArxnidx = find(contains(newrxn.ID,NArxn));
newrxn.GPR(NArxnidx) = [];
newrxn.GPR = strrep(newrxn.GPR,'Saccharomyces_cerevisiae@','');
newrxn.GPRNote(NArxnidx) = [];
newrxn.ID(NArxnidx) = [];
newrxn.Rev(NArxnidx) = [];
newrxn.rxnECNumbers(NArxnidx) = [];
newrxn.rxnKEGGID(NArxnidx) = [];
newrxn.rxnMetaNetXID(NArxnidx) = [];
newrxn.rxnNames(NArxnidx) = [];
newrxn.rxnpathway(NArxnidx) = [];

%check whether one metabolite appear more than once in one rxn
rxncheckcompartment = [];
for i = 1:length(newrxn.ID)
    j = find(strcmp(matrix.rxnIDs,newrxn.ID{i}));
    Met = matrix.metIDs(j);
    [met_new,idx] = unique(Met);
    if length(met_new) < length(Met)
        met_rep = j(setdiff(1:numel(Met), idx));
        if length(met_rep) == 1
            if strcmp(matrix.metIDs(met_rep),'H+')
                met_rep = [];
            end
        end
        if ~isempty(met_rep)
            rxncheckcompartment = [rxncheckcompartment;newrxn.ID(i),met_rep];
        end
    end
end

% change one metablite comp to be extracelluar since it is a transport rxn
for i = 1:length(rxncheckcompartment(:,1))
    j = find(strcmp(matrix.rxnIDs,rxncheckcompartment{i,1}) & strcmp(matrix.mettype,'reactant'));
    Met = matrix.metIDs(j);
    idx = j(ismember(Met,matrix.metIDs(rxncheckcompartment{i,2})));
    matrix.metcompartments(idx) = {'extracellular'};
end
cd otherchanges/
[model,rxnUpdateGPR,EnergyResults,RedoxResults,MassChargeresults] = addPanModelRxn(model,matrix,newmet,newrxn);

%% update the grRules for existing reactions
% We found there are 63 rxns are exisitng in the original model, will check
% that and then decide whthether we should update gpr or not.
model = rmfield(model,'grRules');
model.grRules = printGPRForRxns(model,model.rxns);
[~,ID] = ismember(rxnUpdateGPR(:,1),model.rxns);
oldGPR = printGPRForRxns(model,model.rxns(ID));

% if old GPR is empty, we update the rxn with new GPR and the rxn is essential, it will cause the
% problem of maybe deleting this one will lead to no growth.
% rxnUpdateGPR(cellfun(@isempty,oldGPR),:) = [];
% oldGPR(cellfun(@isempty,oldGPR)) = [];

linker = repmat({' or '},length(oldGPR(:,1)),1);
linker(cellfun(@isempty,oldGPR))= {''};
oldGPR = [oldGPR,linker];
newGPR = [oldGPR,rxnUpdateGPR(:,3)];
newGPR = join(newGPR);
newGPR = strrep(newGPR,'  ',' ');
newGPR = strtrim(newGPR);
mapping.rxnIDs  = rxnUpdateGPR(:,1);
mapping.new_GPR  = newGPR;
model = AddMissingOrthologsInYeastGEM(model,mapping);

%% manual curation for model reversibilities
[model,changes] = ManualCuration(model);
cd ../

%% load matrix information and Generate the specific models
load('ortholog_changedtoPanID.mat') % load all orthologs which will be added back to the specific models
Genes_Query = [model.genes;ortholog(:,1);ortholog(:,2)];
Genes_Query = unique(Genes_Query);

% % update the gpr rules by including orthlogs: A or B to A or B or Panortholog
% model = UpdatePanGPR(ortholog);% cause updating the orthlogs will
% generate a lot of grRules, so we would rather add the ortholog back at
% the specific models section

% load mapping IDs for the variables mappingID.OGIDs are standard IDs in the strain-gene presenceAvsence data,
% mappingID.panIDs are panIDs that we used to construct the panmodel
format = '%s %s';
fID       = fopen('../../find_homolog_for_panID_py/data/representatives.tsv');
mapping  = textscan(fID,format,'Delimiter','\t','HeaderLines',1);
mappingID.OGIDs      = mapping{1};
mappingID.panIDs = mapping{2};
fclose(fID);

%load presenceAvsence data
RecName = '../modelRelated/343_gene_pa_table.csv'; % setting file name
RecStore = datastore(RecName,'ReadVariableNames',true); % set whether should we read variables name
RecStore.Delimiter = ',';

% find columns that match model.genes. There is a problem with the
% variableName, fix the gene name by replaceing '_' to '-' to map the gene
% name in the model
m = RecStore.VariableNames(2:end);
[~,ID] = ismember(m,mappingID.OGIDs);
x = find(ID~=0);
m(x) = mappingID.panIDs(ID(x));
for i = 1:length(m)
    if length(m{i})== 9
        m(i) = strrep(m(i),'_','-');
    end
end

% there is a issue of panID selection, should be fixed using all panID we
% set as in the species specific file
[~,ID] = ismember(geneID_core,m);
m(ID(ID~=0)) = panID_final(ID~=0);
[~,ID] = ismember(Genes_Query,m);
if  ~all(ID)
    warning([Genes_Query{ID ==0},' not in the geneExistence Matrix, please check that.'])
end
selected = find(ID~=0);
selected_mapping = m(ID(selected)); % mapping back to OG name
m =RecStore.VariableNames(2:end);
selected = m(ID(selected)); % matach back the pangenes

RecStore.SelectedVariableNames = selected; % select variables to read
genesMatrix = readall(RecStore);% read selected variables
StrianData.genes = selected_mapping';
StrianData.levels = table2array(genesMatrix)';
RecStore.SelectedVariableNames = RecStore.VariableNames(1);
genesMatrix = readall(RecStore);% read the strain name
StrianData.strains = table2cell(genesMatrix(:,1));
panmodel = model;

clearvars -except model ortholog StrianData panmodel

cd otherchanges/
for i = 1:length(StrianData.strains)
    disp([StrianData.strains{i},' No.',num2str(i)])
    [~,ID] = ismember(StrianData.strains(i),StrianData.strains);
    lvl = StrianData.levels(:,ID);
    lvl_temp = logical(lvl);
    genesExist = StrianData.genes(lvl_temp);
    [~,ID] = ismember(genesExist,ortholog(:,2));
    ortholog_strian = ortholog(ID(ID~=0),:);
    model_temp = model;
    if ~isempty(ortholog_strian)
    model_temp = UpdatePanGPR(ortholog_strian,model);
    end
    [reducedModel,resultfile] = SpecificModel(model_temp,StrianData,StrianData.strains(i),'../../ModelFiles/mat');
end
