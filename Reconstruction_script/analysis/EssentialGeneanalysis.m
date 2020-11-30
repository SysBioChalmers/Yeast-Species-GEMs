function [accurancy,sensitivity,specificity,mcc] = EssentialGeneanalysis(strains,inputpath)
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
exp_result_viable = cell(5000,length(strains));
model_result_viable = cell(5000,length(strains));
exp_result_inviable = cell(5000,length(strains));
model_result_inviable = cell(5000,length(strains));

current_path = pwd;
% load the model
for i = 1:length(strains)
    disp(['No.',num2str(i),'/',num2str(length(strains))])
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

[~,idx] = ismember(exp_viable,upper(model.genes));
exp_viable_prot  = model.proteins(idx);
[~,idx] = ismember(exp_inviable,upper(model.genes));
exp_inviable_prot = model.proteins(idx);
exp_result_viable(1:length(exp_viable_prot),i) = exp_viable_prot;
exp_result_inviable(1:length(exp_inviable_prot),i) = exp_inviable_prot;

%calculate the growth rate after the single gene deletion using the original model and update model
grRatio = singleGeneDeletion(model);
mod_viable  = model.genes(grRatio >= ko_tol);
mod_viable = intersect(upper(mod_viable),verifiedORFs);
mod_inviable = model.genes(grRatio < ko_tol );
mod_inviable = intersect(upper(mod_inviable),verifiedORFs);

mod_viable_prot  = model.proteins(grRatio >= ko_tol);
mod_inviable_prot = model.proteins(grRatio < ko_tol );
model_result_viable(1:length(mod_viable_prot),i) = mod_viable_prot;
model_result_inviable(1:length(mod_inviable_prot),i) = mod_inviable_prot;

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
accurancy(i) = (n_tp+n_tn)/(n_tp+n_tn+n_fn+n_fp)
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


% plot the figure
P = [accurancy;specificity;sensitivity;mcc;n_tps;n_tns;n_fps;n_fns];
writematrix(P,'accuracyrestlt.txt')
save('result_newcomplex.mat')

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
cd ../data/physiology/EssentialGene_ML
% load essential data
fileName = [strain,'.txt'];
fID       = fopen(fileName);
temp  = textscan(fID,'%s%s','Delimiter','\t','HeaderLines',1);
fclose(fID);

%verifiedORFs     = temp{2};
%condition        = temp{3};
verifiedORFs     = temp{1};
condition        = temp{2};
inviableORFs = verifiedORFs(ismember(condition,'E')|ismember(condition,'Essential'));% essential
    verifiedORFs = upper(verifiedORFs);
    verifiedORFs = unique(verifiedORFs);
    inviableORFs = upper(inviableORFs);
    inviableORFs = unique(inviableORFs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
