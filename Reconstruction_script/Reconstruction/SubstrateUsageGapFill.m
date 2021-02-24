function [model_original] = SubstrateUsageGapFill(model_original,strains,inputpath)
% This function is to fill the gaps for the substrateUsage
%inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';
fid2 = fopen('../data/gapfill/substrate_usage_rxn.tsv');
format = repmat('%s ',1,4);
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(temp)
    rxnlist(:,i) = temp{i};
end
%load('StrainData.mat')

% always fill the gap from the panmodel, so that the specific models can
% have the same rxn_ids
cd otherchanges/
[model_original,addedrxn,EnergyResults,MassChargeresults,RedoxResults] = addRxnfromMNX(model_original,rxnlist(2:end,1),'cytoplasm','extracellular');
changes = [];
% manual curation for some reactions
[~,idx] = ismember('MNXR115392',model_original.rxnMetaNetXID); %inulin
rxnformula = 'H2O [extracellular] + inulin [extracellular]  -> 2 D-fructose [extracellular]';
model_original = changerxn(model_original,model_original.rxns{idx},rxnformula);
printRxnFormula(model_original,'rxnAbbrList',model_original.rxns(idx),'metNameFlag',true)
changes = [changes; model_original.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];

[~,idx] = ismember('MNXR95249',model_original.rxnMetaNetXID); %N-Acetyl-D-glucosamine
rxnformula = 'ATP [cytoplasm] + N-Acetyl-D-glucosamine [extracellular]  -> ADP [cytoplasm] + H+ [cytoplasm] + N-acetyl-D-glucosamine&6-phosphate [cytoplasm]';
model_original = changerxn(model_original,model_original.rxns{idx},rxnformula);
printRxnFormula(model_original,'rxnAbbrList',model_original.rxns(idx),'metNameFlag',true)
changes = [changes; model_original.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];

[~,idx] = ismember('MNXR112943',model_original.rxnMetaNetXID); %(R,R)-2,3-butanediol acetoin to R format
rxnformula = 'coenzyme&A [cytoplasm] + NAD [cytoplasm] + (R)-acetoin [cytoplasm]  -> acetaldehyde [cytoplasm] + acetyl-CoA [cytoplasm] + H+ [cytoplasm] + NADH [cytoplasm]';
model_original = changerxn(model_original,model_original.rxns{idx},rxnformula);
printRxnFormula(model_original,'rxnAbbrList',model_original.rxns(idx),'metNameFlag',true)
changes = [changes; model_original.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];

[~,idx] = ismember('MNXR101023',model_original.rxnMetaNetXID); % alpha-D-Glucose to glucose format
rxnformula ='H2O [cytoplasm] + lactose [cytoplasm] -> D-galactose [cytoplasm] + D-glucose [cytoplasm]';
model_original = changerxn(model_original,model_original.rxns{idx},rxnformula);
printRxnFormula(model_original,'rxnAbbrList',model_original.rxns(idx),'metNameFlag',true)
changes = [changes; model_original.rxnMetaNetXID(idx),{'changecomp'},{'for sugar degradation rxns'}];

rxnformula ='H2O [cytoplasm] + beta-cellobiose [cytoplasm] 	->	2 D-glucose [cytoplasm]';
model_original = changerxn(model_original,'r_5112_fwd',rxnformula);

[~,idx1] = ismember('alpha-D-Glucose [cytoplasm]',model_original.metNames);
[~,idx2] = ismember('D-glucose [cytoplasm]',model_original.metNames);
rxnidx = find(model_original.S(idx1,:));
model_original.S(idx2,rxnidx) = model_original.S(idx1,rxnidx);
model_original.S(idx1,rxnidx) = 0;
printRxnFormula(model_original,'rxnAbbrList',model_original.rxns(rxnidx),'metNameFlag',true)

[~,idx1] = ismember('beta-D-Glucose [cytoplasm]',model_original.metNames);
[~,idx2] = ismember('D-glucose [cytoplasm]',model_original.metNames);
rxnidx = find(model_original.S(idx1,:));
model_original.S(idx2,rxnidx) = model_original.S(idx1,rxnidx);
model_original.S(idx1,rxnidx) = 0;
printRxnFormula(model_original,'rxnAbbrList',model_original.rxns(rxnidx),'metNameFlag',true)


[~,idx2] = ismember('D-glucose [cytoplasm]',model_original.metNames);
rxnidx = find(model_original.S(idx2,:));
idx = (sum(full(model_original.S(:,rxnidx))~=0,1) == 1);
model_original = removeReactions(model_original,model_original.rxns(rxnidx(idx)));

%Remove unused metabolites
[usedMets, ~]=find(model_original.S);
unUsedMets=true(numel(model_original.mets),1);
unUsedMets(usedMets)=false;
model_original=removeMets(model_original,unUsedMets,false,false,false,false);

missingMets = {'inulin [extracellular]','beta-cellobiose [extracellular]','Mannitol [extracellular]','D-gluconate [extracellular]','2-dehydro-D-gluconate [extracellular]','ethylamine [extracellular]'};
formulas = {'C12H22O11';'C12H22O11';'C6H14O6';'C6H11O7';'C6H9O7';'C2H7N'};
[~,idx] = ismember(missingMets,model_original.metNames);
model_original.metFormulas(idx(idx~=0)) = formulas(idx~=0);


missingMets = {'N-Acetyl-D-glucosamine [extracellular]','Mannitol [extracellular]','lactose [extracellular]'};
formulas = {'C8H15NO6';'C6H14O6';'C12H22O11'};
[~,idx] = ismember(missingMets,model_original.metNames);
model_original.metFormulas(idx(idx~=0)) = formulas(idx~=0);

% add those rxns for each strain
current_path = pwd;
for i = 2:length(rxnlist(1:end,1))
    strainslist = rxnlist(i,3);
    strainslist = split(strainslist,',');
    strainslist = strrep(strainslist,' ','');
    if strcmp(strainslist,'')
        strainslist = strains;
    end
    strainslist=regexprep(strainslist,'_16....$','');
    [~,idx] = ismember(lower(strainslist),lower(strains));
    [~,rxnIdx] = ismember(rxnlist(i,1),model_original.rxnMetaNetXID);
    for j = 1:length(idx)
        if idx(j) ~=0
            m = strains{idx(j)};
            cd(inputpath)
            reducedModel = load([m,'.mat']);
            reducedModel = reducedModel.reducedModel;
            cd(current_path)
            reducedModel = addrxnBack(reducedModel,model_original,model_original.rxns(rxnIdx),{''});
            cd(inputpath)
            save([m,'.mat'],'reducedModel')
        else
            warning(['no species found for',strainslist{j}])
        end
    end
end

cd(current_path)