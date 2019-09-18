

% This function is the main function for generating panmodel and specific
% models
% Load ortholog information
fid      = fopen('../../find_homolog_for_panID_py/result/pan_hit_mapping_panYeast_v2_PI@70.tsv');
orth     = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
ortholog(:,1)     = orth{1};
ortholog(:,2) = orth{2};
fclose(fid);

% update the gpr rules by including orthlogs: A or B to A or B or Panortholog
model = UpdatePanGPR(ortholog);

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
newrxn.GPR = replace(newrxn.GPR,'NA','');
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
    matrix.metcompartments(idx) = {'extracelluar'};
end
[model,rxnUpdateGPR] = addPanModelRxn(model,matrix,newmet,newrxn);

% We found there are 62 rxns are exisitng in the original model, will check
% that and then decide whthether we should update gpr or not.

