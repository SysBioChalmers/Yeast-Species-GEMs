function accuracy = PhenSpecificRxn(panmodel,rxnMatrix,substrate,aerobic,Subtype)
% calculate one phenotyperelated rxns, compare with the alternative rxns,
% find clade which doesn't inckude this rxn, find percentage of non-used
% species does not contain this rxn.

% This pipeline is to analyze the rxn existence in all species



% load substrate usage data
% find the phenotype
fid2 = fopen('../../data/physiology/Biolog_substrate.tsv');
format = repmat('%s ',1,333);
format = strtrim(format);
substrate = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(substrate)
    data(:,i) = substrate{i};
end
strainlist = data(1,5:end);
data(1,:) = [];
SubBiologNames = data(:,1);
SubModelNames = data(:,2);
SubConditions = data(:,3);
Subtypes = data(:,4);
data(:,1:4) = [];
fclose(fid2);

%load clade information
[~, ~, Strain_information]=xlsread('../data/genome_summary_332_yeasts.xlsx','clades');
Strain_information = Strain_information(2:end,:);
clades = unique(Strain_information(:,2));


newModel = panmodel;
cd ../Reconstruction/otherchanges/
metE = strcat(substrate,' [cytoplasm]');
metC = strcat(substrate,' [cytoplasm]');
[~,metEOrd] = ismember(metE,newModel.metNames);
[~,metCOrd] = ismember(metC,newModel.metNames);
%check whther exchange reaction for this met exists
duplicateE = find(sum(newModel.S ~= 0, 1) == 1 & any(newModel.S == -1, 1) & any(newModel.S(metEOrd(metEOrd~=0), :), 1));
if ~isempty(duplicateE)
    warning(['Model already has the same exchange reaction you tried to add: ', newModel.rxns{duplicateE}]);
    ExchRxn  = newModel.rxns{duplicateE};
    ExchLB   = newModel.lb(duplicateE);
else
    if metEOrd == 0
        newID       = getNewIndex(newModel.mets); % here we make sure that the model
        metE_ID     = strcat('s_',newID,'[e]');
        newModel    = addMetabolite(newModel,metE_ID,'metName',metE);
        [~,metEOrd] = ismember(metE,newModel.metNames);
    end
    %adding exchange rxn for this metE
    newID    = getNewIndex(newModel.rxns);
    ExchRxn  = ['r_' newID];
    ExchLB   = 0;
end
metE_ID = newModel.mets{metEOrd};
[newModel,rxnIDexists] = addReaction(newModel,ExchRxn, ...
    'reactionName', [substrate, ' exchange'], ...
    'metaboliteList', cellstr(metE_ID), 'stoichCoeffList', -1, ...
    'lowerBound', ExchLB, 'upperBound', 1000, 'subSystem', '', ...
    'checkDuplicate', false);
newModel = rmfield(newModel,'grRules');

% Fix confidence score:
SubRxnIndex = findRxnIDs(newModel,ExchRxn);
newModel.rxnConfidenceScores(SubRxnIndex) = NaN;    %exchange rxns
newModel.rxnNotes{SubRxnIndex} = ['NOTES: added after the substrate usage; ', newModel.rxnNames{SubRxnIndex}];

% Change media and test:
newModel_test = newModel;
newModel_test = minimal_Y6(newModel_test);
glcExchange = {'r_1714'}; % D-glucose exchange
amoExchange = {'r_1654'}; % ammonium exchange
phoExchange = {'r_2005'}; % phosphate exchange
sulExchange = {'r_2060'}; % phosphate exchange
amoExchangeIndex = findRxnIDs(newModel_test,amoExchange);
phoExchangeIndex = findRxnIDs(newModel_test,phoExchange);
glcExchangeIndex = findRxnIDs(newModel_test,glcExchange);
sulExchangeIndex = findRxnIDs(newModel_test,sulExchange);
if Subtype == 'C'
    newModel_test.lb(glcExchangeIndex) = 0;
elseif Subtype == 'N'
    newModel_test.lb(amoExchangeIndex) = 0;
elseif Subtype == 'P'
    newModel_test.lb(phoExchangeIndex) = 0;
elseif Subtype == 'S'
    newModel_test.lb(sulExchangeIndex) = 0;
end
newModel_test.lb(SubRxnIndex) = -10;
if strcmp(aerobic,'0') || strcmp(aerobic,'anaerobic')% which means fermentaion
    newModel_test = anaerobicModel(newModel_test);
end

% Simulate model: glpk has a problem, but cobra toolbox doesn't support
% gurobi as sovler for matlab2019b, so we will use RAVEN function here, but
% we need the field model.rev
for i=1:numel(newModel_test.rxns)
    if newModel_test.lb(i)<0
        newModel_test.rev(i,1)=1;
    else
        newModel_test.rev(i,1)=0;
    end
end
for i=1:numel(model.rxns)
    if model.lb(i)<0
        model.rev(i,1)=1;
    else
        model.rev(i,1)=0;
    end
end

sol = solveLP(newModel_test,1);
sol_glc = solveLP(model,1);

clearvars -except model panmodel substrate Subtype aerobic sol sol_glc Strain_information cladesdata Subtypes SubModelNames SubConditions rxnMatrix

%sol = optimizeCbModel(newModel_test,'max','one');
%sol_glc = optimizeCbModel(panmodel,'max','one');

rxn_specfic = setdiff(newModel_test.rxns(sol.x~=0),newModel_test.rxns(sol_glc.x~=0));

% find the rxn that are in the group
cd ../
alterrxns = model.rxns(sum(rxnMatrix,1) ~= 343);
rxn_specfic = intersect(rxn_specfic,alterrxns);

sub_strain_model = strains((rxnMatrix(:,783) ~= 0))


% find the species experimentally 
idx = find(ismember(SubModelNames,substrate) & ismember(Subtypes,subtype) & ismember(SubConditions,aerobic)) ;
sub_strain_exp = strainlist(strcmp(data(idx,:),'1')| strcmp(data(idx,:),'v'));

  



