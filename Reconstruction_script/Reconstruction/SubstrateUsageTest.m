function [accuracy] = SubstrateUsage(strains,model_original,panmodelTest,inputpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubstrateUsage
%
% Automatically adds exchange reactions for every metabolite in
% ComplementaryData/physiology/Biolog_substrate.tsv and checks whether
% it can be used as a "solo" substrate.
%
% NOTE: requires COBRA
%
% Feiran Li     2018-08-25
% Feiran Li     2019-12-12 -Modify that to fit the Strain specific models;
% change that into a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load data:
fid2 = fopen('../data/physiology/Biolog_substrate.tsv');
format = repmat('%s ',1,333);
format = strtrim(format);
substrate = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(substrate)
    data(:,i) = substrate{i};
end
strainlist = data(1,5:end);
data(1,:) = [];
SubBiologName = data(:,1);
SubModelName = data(:,2);
SubCondition = data(:,3);
Subtype = data(:,4);
data(:,1:4) = [];
fclose(fid2);

path = pwd;



%% test for panmodel to see wether the panmodel can ultilize those metabolites or not.
% Main loop:
cd otherChanges
if strcmp(panmodelTest,'true') || panmodelTest
    ExchRxn             = '';
    TransRxn            = '';
    GapfillMets         = '';
    TransEotherfillMets = '';
    TransECfillMets     = '';
    FBAresult           = '';
    for i = 1:length(SubBiologName)
        newModel = model_original;
        if ~isequal(SubModelName{i},'')
            metE = strcat(SubModelName{i},' [extracellular]');
            metC = strcat(SubModelName{i},' [cytoplasm]');
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
                'reactionName', [SubModelName{i}, ' exchange'], ...
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
            if Subtype{i} == 'C'
                newModel_test.lb(glcExchangeIndex) = 0;
            elseif Subtype{i} == 'N'
                newModel_test.lb(amoExchangeIndex) = 0;
            elseif Subtype{i} == 'P'
                newModel_test.lb(phoExchangeIndex) = 0;
            elseif Subtype{i} == 'S'
                newModel_test.lb(sulExchangeIndex) = 0;
            end
            newModel_test.lb(SubRxnIndex) = -10;
            if strcmp(SubCondition{i},'0') % which means fermentaion
                newModel_test = setParam(newModel_test,'obj','r_1761', 1);
                newModel_test = setParam(newModel_test,'eq','r_2111', 0.01);
            end
            
            % Simulate model:
            sol = optimizeCbModel(newModel_test);
            if sol.obj > 0.000001
                model_original = newModel;
                fprintf(metE,'can be used as solo substrate\n');
                FBAresult = [FBAresult;{'1'}];
            else
                FBAresult = [FBAresult;{'0'}];
            end
        else
            FBAresult = [FBAresult;{'n'}];
        end
    end
end
%% Test the ssModels subtrate usage information
for k = 1:length(strainlist)
    cd(inputpath)
    load([strainlist{k},'.mat'])
    model = reducedModel;
    cd(current_path)
    for i = 1:length(SubBiologName)
        newModel = model;
        if ~isequal(SubModelName{i},'')
            metE = strcat(SubModelName{i},' [extracellular]');
            metC = strcat(SubModelName{i},' [cytoplasm]');
            [~,metEOrd] = ismember(metE,newModel.metNames);
            %check whther exchange reaction for this met exists
            duplicateE = find(sum(newModel.S ~= 0, 1) == 1 & any(newModel.S == -1, 1) & any(newModel.S(metEOrd(metEOrd~=0), :), 1));
            if ~isempty(duplicateE)
                warning(['Model already has the same exchange reaction you tried to add: ', newModel.rxns{duplicateE}]);
                ExchRxn  = newModel.rxns{duplicateE};
                SubRxnIndex = findRxnIDs(newModel,ExchRxn); % find the Exchange rxn for that metE.
            end
            % add excahnge reactions
            [~,rxnIdx] = ismember([SubModelName{i}, ' exchange'],model_original.rxnNames);
            if rxnIdx ~= 0 
            newModel = addrxnBack(newModel,model_original,model_original.rxns(rxnIdx),{''});
            SubRxnIndex = findRxnIDs(newModel,model_original.rxns(rxnIdx)); % find the Exchange rxn for that metE.

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
            if Subtype{i} == 'C'
                newModel_test.lb(glcExchangeIndex) = 0;
            elseif Subtype{i} == 'N'
                newModel_test.lb(amoExchangeIndex) = 0;
            elseif Subtype{i} == 'P'
                newModel_test.lb(phoExchangeIndex) = 0;
            elseif Subtype{i} == 'S'
                newModel_test.lb(sulExchangeIndex) = 0;
            end
            newModel_test.lb(SubRxnIndex) = -10;
            if strcmp(SubCondition,'0') % which means fermentaion
                newModel_test = setParam(newModel_test,'obj','r_1761', 1);
                newModel_test = setParam(newModel_test,'eq','r_2111', 0.01);
            end
            
            % Simulate model:
            sol = optimizeCbModel(newModel_test);
            if sol.obj > 0.000001
                model = newModel;
                fprintf(metE,'can be used as solo substrate\n');
                FBAresult(i,k) = {'1'};
            else
                FBAresult(i,k) = {'0'};
            end
            end
        else
            FBAresult(i,k) = {'n'};
        end
    
    reducedModel = model;
    %cd(outputpath)
    save([strainlist{k},'.mat'],'reducedModel')
    end
end


% plot the difference
for i = 1:length(strainlist)
    nocheck = find(strcmp(FBAresult(:,i),'n'));
    tp(i) = length(intersect(find(strcmp(FBAresult(:,i),'1')),setdiff(find(strcmp(data(:,i),'1')|strcmp(data(:,i),'v')),nocheck)));
    tn(i) = length(intersect(find(strcmp(FBAresult(:,i),'0')),setdiff(find(strcmp(data(:,i),'0')),nocheck)));
    fp(i) = length(intersect(find(strcmp(FBAresult(:,i),'1')),setdiff(find(strcmp(data(:,i),'0')),nocheck)));
    fn(i) = length(intersect(find(strcmp(FBAresult(:,i),'0')),setdiff(find(strcmp(data(:,i),'1')|strcmp(data(:,i),'v')),nocheck)));
    accuracy(i) = (tp(i) + tn(i))/(tp(i) + tn(i) + fn(i) + fp(i));
end
% using clade information to collect this
h = cdfplot(accuracy);
set( h, 'linewidth',3);
set(gca,'FontSize',20,'FontName','Helvetica');

ylabel('Percentage','FontSize',24,'FontName','Helvetica','Color','k');
xlabel('Substrate prediction accuracy','FontSize',24,'FontName','Helvetica','Color','k');

% based on clade
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
Strain_information(:,i) = data{i};
end
fclose(fid2);
clades = unique(Strain_information(:,2));
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};

clades(1) = [];
group = [];
clade_av = [];
for i = 1:length(clades)
idx = ismember(Strain_information(:,2),clades(i));
[~,ID] = ismember(Strain_information(idx,1),strainlist);
clade_av = [clade_av;accuracy(ID(ID~=0))'];
group = [group;(i-1)*ones(length(accuracy(ID(ID~=0))),1)];
end
h = boxplot(clade_av,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Substrate prediction ccuracy','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',60)
set(gcf,'position',[200 0 350 190]);
set(gca,'position',[0.11 0.3 0.87 0.75]);


% plot the difference
for i = 1:75
    nocheck = find(strcmp(FBAresult(i,:),'n'));
    tp(i) = length(intersect(find(strcmp(FBAresult(i,:),'1')),setdiff(find(strcmp(data(i,:),'1')|strcmp(data(i,:),'v')),nocheck)));
    tn(i) = length(intersect(find(strcmp(FBAresult(i,:),'0')),setdiff(find(strcmp(data(i,:),'0')),nocheck)));
    fp(i) = length(intersect(find(strcmp(FBAresult(i,:),'1')),setdiff(find(strcmp(data(i,:),'0')),nocheck)));
    fn(i) = length(intersect(find(strcmp(FBAresult(i,:),'0')),setdiff(find(strcmp(data(i,:),'1')|strcmp(data(i,:),'v')),nocheck)));
    acc_sub(i) = (tp(i) + tn(i))/(tp(i) + tn(i) + fn(i) + fp(i));
end
hold on
sub = SubModelName;
sub(isnan(acc_sub)) = [];
acc_sub(isnan(acc_sub)) = [];
bar([1:1:length(sub)],acc_sub)
set(gca,'XTick',1:1:length(sub));
set(gca,'XTickLabel',sub);
xtickangle(90);
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Accuracy','FontSize',20,'FontName','Helvetica','Color','k');
xlabel('Substrate','FontSize',20,'FontName','Helvetica','Color','k');
%% fix the gap between the panmodel and the ssModels.


Update results from model:
cd ..
fid2 = fopen('../ComplementaryData/physiology/Biolog_substrate.tsv','w');
formatSpec = '%s\t%s\t%s\t%s\t%s\n';
fprintf(fid2,formatSpec,'Substrate','Name_in_Model','Substrate_type','Growth_Biolog','Growth_Model');
for i = 1:length(FBAresult)
    fprintf(fid2,formatSpec,char(FBAresult(i,1)),char(FBAresult(i,2)),char(FBAresult(i,3)),char(FBAresult(i,4)),char(FBAresult(i,5)));
end

% Save model:
saveYeastModel(model)
cd modelCuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
