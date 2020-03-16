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

% the direction cannot uptake exthanol
[~,idx] = ismember('MNXR101464',model.rxnMetaNetXID);
model.lb(idx) = -1000;
printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true)

changes = [changes; model.rxnMetaNetXID(idx),{'changerev'},{'for uptake of methanol'}];
