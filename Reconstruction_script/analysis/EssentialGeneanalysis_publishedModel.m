function [accurancy,sensitivity,specificity,mcc] = EssentialGeneanalysis_publishedModel(strains,inputpath)
% This function is to perform essential analysis and plot the figure

% essential data is loaded from the data filefolder

% essentialGenes
%   Modify media + find essential genes in model. Adapted from:
%   https://doi.org/10.1371/journal.pcbi.1004530
%
%   accurancy   (tp+tn)/(tp+tn+fn+fp)
%   tp          true positives
%   tn          true negatives
%   fn          false negatives
%   fp          false positives
%
%
%   Feiran Li, 2019-09-25
%   Feiran Li, 2020-04-08
%

%initCobraToolbox

ko_tol = 1e-6;

current_path = pwd;
% load the model
for i = 1:length(strains)
    strain = strains{i};
    cd(inputpath)
    model = load([strain,'.mat']);
    model = model.reducedModel;
%constraints from genotype: check the genotype of the strains used in deletion experiment
    % no model changes based on genotype
cd(current_path)
%constraints from growth medium:
model = complete_Y7(model); %change Y7 model medium to the Kennedy synthetic complete medium

%import the gene deletion resuts:
    % compare gene essentiality predictions in minimal media with glucose
    % as sole carbon source against a list of "essential" genes. The
    % problem with this is that the "essential" genes are defined from the
    % stanford yeast knock-out collection, which was screened in complex
    % media supplemented with amino acids required by genetic markers, so
    % it's not an ideal reference gene list. However, we can compare all
    % models against it, so it's useful for comparative purposes. The
    % reference lists of genes are at the end of this function.
[inviableORFsAll,verifiedORFs] = loaddata(strain);
exp_inviable = intersect(upper(model.genes),inviableORFsAll);
exp_inviable = intersect(exp_inviable,verifiedORFs);
exp_viable = setdiff(upper(model.genes),inviableORFsAll);
exp_viable = intersect(exp_viable,verifiedORFs);

%calculate the growth rate after the single gene deletion using the original model and update model
grRatio = singleGeneDeletion(model);
mod_viable  = model.genes(grRatio >= ko_tol);
mod_viable = intersect(upper(mod_viable),verifiedORFs);
mod_inviable = model.genes(grRatio < ko_tol );
mod_inviable = intersect(upper(mod_inviable),verifiedORFs);

%obtain the essential gene and non-essential gene based on the above two calculations
%calculate the prediction number in the followed four type by comparing the
%prediction and experiment
%TP: true positive; TN: true negative; FN: false negative; FP: false positive
tp = intersect(exp_viable,mod_viable); n_tp = length(tp);
tn = intersect(exp_inviable,mod_inviable); n_tn = length(tn);
fp = intersect(exp_inviable,mod_viable); n_fp = length(fp);
fn = intersect(exp_viable,mod_inviable); n_fn = length(fn);

%compare the prediction performances of two models
%prediction accuracy was used to evaluate the quality of model update in
%each step
accurancy(i,1) = (n_tp+n_tn)/(n_tp+n_tn+n_fn+n_fp)
sensitivity(i) = (n_tp/(n_tp+n_fn));
specificity(i) = (n_tn/(n_tn+n_fp));
positivePredictive(i) = (n_tp/(n_tp+n_fp));
negativePredictive(i) = (n_tn/(n_fn+n_tn));
mcc(i) = (n_tp * n_tn - n_fp * n_fn)/...
sqrt((n_tp + n_fp)*(n_tp + n_fn)*(n_tn + n_fp)*(n_tn + n_fn));
%geoMean = (sensitivity * specificity)^.5;

n_tps(i) = n_tp;
n_tns(i) = n_tn;
n_fps(i) = n_fp;
n_fns(i) = n_fn;
acc(i) = accurancy(i);
end

% combine with published model result
% load mapping list published model/species
fid2 = fopen('../data/PublishedModel_list.tsv');
format = '%s %s %s %s %s %s %s';
tmp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(tmp)
data(:,i) = tmp{i};
end
fclose(fid2);
[~,species_idx] = ismember(strains,data(:,2));
accuray_published = zeros(length(strains),1);
accuray_published(species_idx~=0) = cellfun(@str2num,data(species_idx(species_idx~=0),6));
accuracy_final = [accurancy,accuray_published];
h = bar(accuracy_final);
h(2).FaceColor = [178,24,43]/255;
h(1).FaceColor = [33,102,172]/255;
xticklabels(strrep(strains,'_',' '))
xtickangle(45);
legend({'Models this work','Published models'})
yticks([0:0.2:1])
ylim([0,1])
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Essential gene prediction accuracy','FontSize',24,'FontName','Helvetica','Color','k');

% plot the figure
P = [accurancy;specificity;sensitivity;mcc];
spider_plot_R2019b(P',...
'AxesInterval', 4,...
'AxesDisplay', 'one',...
'AxesPrecision', 2,...
'AxesLimits', [0, 0, 0, 0; 1, 1, 1, 1],...
'FillOption', 'on',...
'FillTransparency', 0.1,...
'LineWidth', 4,...
'Marker', 'none',...
'AxesLabels', {'accuracy', 'specificity', 'sensitivity', 'mcc'},...
'AxesFontSize', 14,...
'LabelFontSize', 10);

legend(replace(strains,'_',' '), 'Location', 'southoutside');
set(gcf,'position',[200 400 300 200]);



function model = complete_Y7(model)
    % change Y7 model medium to the Kennedy synthetic complete medium
    % NOTE: Y7 doesn't have exchange reactions for dihydrofolate,
    % octadecenoate, octadecynoate, or linoleic acid.

    % start with a clean slate: set all exchange reactions to upper bound =
    % 1000 and lower bound = 0 (ie, unconstrained excretion, no uptake)
    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;

    constrainedUptake = {'r_1604';'r_1639';'r_1873';'r_1879';'r_1880';...
        'r_1881';'r_1671';'r_1883';'r_1757';'r_1891';'r_1889';'r_1810';...
        'r_1993';'r_1893';'r_1897';'r_1947';'r_1899';'r_1900';'r_1902';...
        'r_1967';...
        %'r_1918';... Y7 doesn't have a linoleic acid exchange rxn, which
        %was substituted for octadecenoate and octadecynoate for Y5 and Y6,
        %so leave out.
        'r_1903';'r_1548';'r_1904';'r_2028';...
        'r_2038';'r_1906';'r_2067';'r_1911';'r_1912';'r_1913';'r_2090';...
        'r_1914';'r_2106'};   % potassium exchange};
    glucoseExchange = {'r_1714';};
    unconstrainedUptake = {'r_1672';'r_1654'; ... % ammonium exchange
                    'r_1992'; ... % oxygen exchange
                    'r_2005'; ... % phosphate exchange
                    'r_2060'; ... % sulphate exchange
                    'r_1861'; ... % iron exchange, for test of expanded biomass def
                    'r_1832'; ... % hydrogen exchange
                    'r_2100'; ... % water exchange
                    'r_4593'; ... % chloride exchange
                    'r_4595'; ... % Mn(2+) exchange
                    'r_4596'; ... % Zn(2+) exchange
                    'r_4597'; ... % Mg(2+) exchange
                    'r_2049'; ... % sodium exchange
                    'r_4594'; ... % Cu(2+) exchange
                    'r_4600'; ... % Ca(2+) exchange
                    'r_2020' };

    constrainedUptakeRxnIndexes = findRxnIDs(model,constrainedUptake);
    glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
    unconstrainedUptakeRxnIndexes = findRxnIDs(model,unconstrainedUptake);
    model.lb(constrainedUptakeRxnIndexes) = -0.5;
    model.lb(glucoseExchangeIndex) = -20;
    model.lb(unconstrainedUptakeRxnIndexes) = -1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inviableORFs,verifiedORFs] = loaddata(strain)
cd ../data/physiology/EssentialGene/
% load essential data
fileName = [strain,'.csv'];
fID       = fopen(fileName);
temp  = textscan(fID,'%s%s%s','Delimiter','\t','HeaderLines',1);
fclose(fID);

verifiedORFs     = temp{2}; 
condition        = temp{3};
inviableORFs = verifiedORFs(ismember(condition,'E'));% essential
    verifiedORFs = upper(verifiedORFs);
    verifiedORFs = unique(verifiedORFs);
    inviableORFs = upper(inviableORFs);
    inviableORFs = unique(inviableORFs);    
    cd ../../../analysis/
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end