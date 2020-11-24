function [model,changes] = ManualCuration(model)
% the model can produce infinite ATP and NADHs
% input is panmodel
changes = [];
% This one is to curate the model

% according to the RXN-5603 in metacyc rxn, so we modify the reversibility
% of this rxn to be irresversible and change the rxn compartment to be in
% mitochondria 
[~,idx] = ismember('MNXR107662',model.rxnMetaNetXID);
model.lb(idx) = 0;
model.S(:,idx) = -model.S(:,idx);
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'RXN-5603'}];
%  according to rxn RXN-14138 in metacyc rxn
[~,idx] = ismember('MNXR112869',model.rxnMetaNetXID);
model.lb(idx) = 0;
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},'RXN-14138'];

% according to rxn r_0658 to curate r_0659 reversibility

[~,idx] = ismember('r_0659',model.rxns);
model.lb(idx) = 0;
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'r_0658'}];


% according to rxn FORMYLTHFDEFORMYL-RXN in metacyc to curate reversibility

[~,idx] = ismember('MNXR99669',model.rxnMetaNetXID);
model.lb(idx) = 0;
model.S(:,idx) = -model.S(:,idx);
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'FORMYLTHFDEFORMYL-RXN'}];

% according to rxn r_0201 in model to curate reversibility

[~,idx] = ismember('MNXR94991',model.rxnMetaNetXID);
model.lb(idx) = 0;
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'r_0201'}];

% There is no direction in metacyc and seed but this one will cause the
% loop of ATP production
[~,idx] = ismember('MNXR106342',model.rxnMetaNetXID);
model.lb(idx) = 0;
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'reason not found'}];

% the direction cannot uptake ,methanol
[~,idx] = ismember('MNXR101464',model.rxnMetaNetXID);
model.lb(idx) = -1000;
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'for uptake of methanol'}];

% change rxn to match glucose in the model
[~,idx] = ismember('MNXR104287',model.rxnMetaNetXID); % glucitol degradate to glucose
rxnformula = 'D-glucitol [cytoplasm] + NADP(+) [cytoplasm]  <=> H+ [cytoplasm] + NADPH [cytoplasm] + D-Glucose [cytoplasm]';
model.rxnMetaNetXID(idx) = {'MNXR104286'};
model = changerxn(model,model.rxns{idx},rxnformula);
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changemet'},{'for using the same glucose'}];

% change rxn compartment
[~,idx] = ismember('MNXR97317',model.rxnMetaNetXID); % melibiose
rxnformula = 'H2O [extracellular] + melibiose [extracellular]  => D-galactose [extracellular] + D-glucose [extracellular]';
model = changerxn(model,model.rxns{idx},rxnformula);
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];

% change met name
[~,idx] = ismember('MNXR104286',model.rxnMetaNetXID); % glucitol
rxnformula = 'D-glucitol [cytoplasm] + NADP(+) [cytoplasm] 	<=>	H+ [cytoplasm] + NADPH [cytoplasm] + D-glucose [cytoplasm]'; 
model = changerxn(model,model.rxns{idx},rxnformula);
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];

% change met name
[~,idx] = ismember('MNXR104286',model.rxnMetaNetXID); % glucitol
rxnformula = 'D-glucitol [cytoplasm] + NADP(+) [cytoplasm] 	<=>	H+ [cytoplasm] + NADPH [cytoplasm] + D-glucose [cytoplasm]'; 
model = changerxn(model,model.rxns{idx},rxnformula);
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];

[~,idx] = ismember('MNXR96241',model.rxnMetaNetXID); % beta-cellobiose
oldgpr = model.grRules(idx);
model = removeReactions(model,model.rxns(idx),'removeUnusedMets',true);
[~,idx] = ismember('MNXR106343',model.rxnMetaNetXID); % beta-cellobiose
newgpr = join([oldgpr,' or ',model.grRules(idx)]);
rxnformula = 'H2O [cytoplasm] + beta-cellobiose [cytoplasm]  <=> 2 D-glucose [cytoplasm]'; 
model = changerxn(model,model.rxns{idx},rxnformula,cell2mat(newgpr));
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)

[~,idx] = ismember('MNXR106516',model.rxnMetaNetXID); % ATP generation cycle
rxnformula = 'fumarate [mitochondrion] + H+ [mitochondrion] + NADH [mitochondrion] => NAD [mitochondrion] + succinate [mitochondrion]'; 
model = changerxn(model,model.rxns{idx},rxnformula);
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)

% curate the grRule for rxn r_4042 raffinose degradation
[~,idx] = ismember('MNXR103420',model.rxnMetaNetXID);
rxnformula = 'raffinose [extracellular] -> D-fructose [extracellular] + melibiose [extracellular]';
model = changerxn(model,model.rxns{idx},rxnformula,'YIL162W');
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'r_4042'}];

% curate the grRule for rxn r_4491 for D_arabionse degradation the same as
% xylose XR(PMID: 27400037)
[~,idx] = ismember('r_4491',model.rxns);
rxnformula = 'D-arabinose [cytoplasm] + H+ [cytoplasm] + NADPH [cytoplasm] => NADP(+) [cytoplasm] + D-arabinitol [cytoplasm]';
model = changerxn(model,model.rxns{idx},rxnformula,'YHR104W or yHAB160_Kazachstania_kunashirensis@Seq_5180 or yHMPu5000026142_Citeromyces_matritensis@Seq_4894 or yHMPu5000035046_Barnettozyma_populi@Seq_3049');
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'r_4491'}];

% get gene name
fid = fopen('../../data/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

for i = 1: length(model.genes)
geneIndex = strcmp(yeast_gene_annotation{1}, model.genes{i});
if sum(geneIndex) == 1 && ~isempty(yeast_gene_annotation{2}{geneIndex})
model.geneNames{i} = yeast_gene_annotation{2}{geneIndex};
else
model.geneNames{i} = model.genes{i};
end
end
for i = 1:length(model.genes)
model.proteins{i} = strcat('COBRAProtein',num2str(i));
end
model = removeUnusedGenes(model);

[~,idx] = ismember('r_0465',model.rxns);
model.rxnKEGGID(idx) = {'R00765'}; % for later gapfilling usage using draftmodels

% [~,idx] = ismember('MNXR101023',model.rxnMetaNetXID); %lactose
% rxnformula = 'H2O [cytoplasm] + lactose [cytoplasm]  <=> D-galactose [cytoplasm] + D-glucose [cytoplasm]';
% model = changerxn(model,model.rxns{idx},rxnformula);
% printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
% changes = [changes; model.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];
% 
% 
% [~,idx] = ismember('MNXR115736',model.rxnMetaNetXID); %lactose
% rxnformula ='H2O [cytoplasm] + lactose [cytoplasm]  => D-glucose [cytoplasm] + beta-D-galactose [cytoplasm]';
% model = changerxn(model,model.rxns{idx},rxnformula);
% printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
% changes = [changes; model.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];
% 
% [~,idx] = ismember('MNXR101000',model.rxnMetaNetXID); %lactose
% rxnformula = 'D-galactose [cytoplasm] + D-glucose [cytoplasm]  => H2O [cytoplasm] + lactose [cytoplasm]';
% model = changerxn(model,model.rxns{idx},rxnformula);
% printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)
% changes = [changes; model.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];
% 

