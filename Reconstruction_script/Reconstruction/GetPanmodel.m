% This function is the main function for generating panmodel and specific
% models



% load new rxn and new metabolites and generate three tsv files for next step: adding new rxns and mets into the model
% mapping metaNetIDs 
format = '%s %s %s %s %s %s %s %s %s %s ';
fID       = fopen('../../Reaction_and_metabolite_standardization/result/new_met_information_from_MNX_database.txt');
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

fID       = fopen('../../Reaction_and_metabolite_standardization/result/new_rxn_information_from_MNX_database.txt');
rxnData = textscan(fID,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newrxn.ID  = rxnData{1};
newrxn.Rev = cellfun(@str2num,replace(rxnData{24},'need_manual_check','1'));
newrxn.GPR = replace(rxnData{8},';',' or ');
newrxn.GPR = replace(newrxn.GPR,'NA',',');
newrxn.GPR = replace(newrxn.GPR,'"','');
newrxn.rxnNames     = rxnData{12};
newrxn.rxnNamesKEGG     = rxnData{11};
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
newrxn.ID(NArxnidx) = [];
newrxn.Rev(NArxnidx) = [];
newrxn.rxnECNumbers(NArxnidx) = [];
newrxn.rxnKEGGID(NArxnidx) = [];
newrxn.rxnMetaNetXID(NArxnidx) = [];
newrxn.rxnNames(NArxnidx) = [];

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
[model,rxnUpdateGPR] = addPanModelRxn(model,matrix,newmet,newrxn);
% We found there are 62 rxns are exisitng in the original model, will check
% that and then decide whthether we should update gpr or not.
model = rmfield(model,'grRules');
model.grRules = printGPRForRxns(model,model.rxns);
[~,ID] = ismember(rxnUpdateGPR(:,1),model.rxns);
oldGPR = printGPRForRxns(model,model.rxns(ID));
% if old GPR is empty, we update the rxn with new GPR and the rxn is essential, it will cause the
% problem of maybe deleting this one will lead to no growth.
rxnUpdateGPR(cellfun(@isempty,oldGPR),:) = [];
oldGPR(cellfun(@isempty,oldGPR)) = [];

linker = repmat({' or '},length(oldGPR(:,1)),1);
linker(cellfun(@isempty,oldGPR))= {''};
oldGPR = [oldGPR,linker];
newGPR = [oldGPR,rxnUpdateGPR(:,3)];
newGPR = join(newGPR);
newGPR = strtrim(newGPR);
mapping.rxnIDs  = rxnUpdateGPR(:,1);
mapping.new_GPR  = newGPR;
model = AddMissingOrthologsInYeastGEM(model,mapping);




% Load ortholog information
fid      = fopen('../../find_homolog_for_panID_py/result/pan_hit_mapping_panYeast_v2_PI@70.tsv');
orth     = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
ortholog(:,1)     = orth{1};
ortholog(:,2) = orth{2};
fclose(fid);

Genes_Query = [model.genes;ortholog(:,1);ortholog(:,2)];
Genes_Query = unique(Genes_Query);
% % update the gpr rules by including orthlogs: A or B to A or B or Panortholog
% model = UpdatePanGPR(ortholog);

% load mapping IDs for the variables mappingID.OGIDs are standard IDs in the strain-gene presenceAvsence data,
% mappingID.panIDs are panIDs that we used to construct the panmodel
format = '%s %s';
fID       = fopen('../../find_homolog_for_panID_py/data/representatives.tsv');
mapping  = textscan(fID,format,'Delimiter','\t','HeaderLines',1);
mappingID.OGIDs      = mapping{1};
mappingID.panIDs = mapping{2};
fclose(fID);

%load presenceAvsence data
RecName = '../gene_pa_table.csv'; % setting file name
RecStore = datastore(RecName,'ReadVariableNames',true); % set whether should we read variables name
RecStore.Delimiter = ','; 

% find columns that match model.genes. There is a problem with the
% variableN 
m = RecStore.VariableNames(2:end);
[~,ID] = ismember(m,mappingID.OGIDs);
x = find(ID~=0);
m(x) = mappingID.panIDs(ID(x));
for i = 1:length(m)
    if length(m{i})== 9
        m(i) = strrep(m(i),'_','-');
    end
end
[~,ID] = ismember(Genes_Query,m);
selected = find(ID~=0);
selected_mapping = m(ID(selected));
m =RecStore.VariableNames(2:end);
selected = m(ID(selected));

RecStore.SelectedVariableNames = selected; % select variables to read
genesMatrix = readall(RecStore);% read selected variables
StrianData.genes = selected_mapping';
StrianData.levels = table2array(genesMatrix)';
RecStore.SelectedVariableNames = RecStore.VariableNames(1);
genesMatrix = readall(RecStore);% read the strain name
StrianData.strains = table2cell(genesMatrix(:,1));
rxnMatrix = zeros(length(model.rxns),1);
panmodel = model;
for i = 1:length(StrianData.strains)
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
    [reducedModel,resultfile] = SpecificModel(model_temp,StrianData,StrianData.strains(i));
    [index,~] = ismember(model.rxns,reducedModel.rxns);
    rxnMatrix(i,:) = transpose(index);
end
    
% analyse the mets that cannot be produced

proMarix = [];
for i = 1:length(StrianData.strains)
    m = StrianData.strains{i};
    load([m,'.mat'])
    model = reducedModel;
    
    bio_rxn = {'r_4048';'r_4049';'r_4050';'r_4598';'r_4599'};% all biomass pseudoreactions except protein, we will manually add amino acid production, since the precursor in protein_pseudoreaction are aa_chargerd tRNAs
    [~,bio_rxn_index] = ismember(bio_rxn,model.rxns);
    mets_test = [];
    for j = 1:length(bio_rxn_index)
        mets_temp = find(model.S(:,bio_rxn_index(j))< 0 );
        mets_test = [mets_test;model.mets(mets_temp)];
        
    end
    
    aa = {'s_1267','s_0956','s_0966','s_0970','s_0974','s_0982','s_1000','s_0992','s_1004','s_1007','s_1017','s_1022','s_1026','s_1030','s_1033','s_1036','s_1041','s_1046','s_1049','s_1052','s_1057'};
%     [~,aa_index] = ismember(aa,model.mets);
%     aa_index = transpose(aa_index);
    mets_test = [mets_test;transpose(aa)];
    
    sol = solveLP(model);
    if ~isempty(sol.x)
    results = canProduce(model,mets_test);
    proMarix(i,:) = transpose(results);
    end
end

% find rxn related to production of that met % ranMatrix stands for rxn
% existence information in each strain, all stores information whether the
% products can be produced or not
% we sort them into three groups: can produce this met and other mets
% tested
mapping = [];
for i = 1:length(proMarix(1,:))
    i
    % group 1 can produce this met
    group1 = StrianData.strains(find(proMarix(:,i))); 
    group1_comrxn = rxnMatrix(find(proMarix(:,i)),:);
    group1_comrxn = panmodel.rxns(all(group1_comrxn,1));
    
    % group 2 cannot produce this met & 90% other mets tested can be
    % produced
    group2 = StrianData.strains(find(~proMarix(:,i))); 
    group2_comrxn = rxnMatrix(find(~proMarix(:,i)),:);
    group2_comrxn = panmodel.rxns(all(group2_comrxn,1));
    
    rxn_test = setdiff(group1_comrxn,group2_comrxn);

    model_out = panmodel;
    for j = 1:length(rxn_test)
        model = model_out;
        model = removeRxns(model, rxn_test{j});
        [model1, rxns]=addExchangeRxns(model,'out',mets_test(i));
        model1 = changeObjective(model1,model1.rxns(end),+1);
        sol = optimizeCbModel(model1,'max','one');
        %sol.f
        if sol.f <= 0
            [~,rxn_index] = ismember(rxn_test(j),panmodel.rxns);
            mapping = [mapping;rxn_test(j),mets_test(i),{strjoin(StrianData.strains(~rxnMatrix(:,rxn_index)),',')},{strcat(num2str(proMarix(~rxnMatrix(:,rxn_index),i)))}];
        else
            model_out = model;
        end
    end
end
[~,rxn_index] = ismember(mapping(:,1),panmodel.rxns);
% Change back to cobra format, in order to get ECnumber,MNXids.
model1 = ravenCobraWrapper(panmodel);
mapping(:,5) = model1.rxnECNumbers(rxn_index);
mapping(:,6) = model1.rxnMetaNetXID(rxn_index);
% This step could take a long time
model1 = UpdatePanGPR(ortholog,panmodel);
mapping(:,7) = model1.grRules(rxn_index);

clearvars -EXCEPT mapping StrianData panmodel
%% This step is to find whether those rxns we found in last step exist in draft model or not


