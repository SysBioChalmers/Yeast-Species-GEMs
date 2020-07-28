% This function is to fill the gaps for the substrateUsage
inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';
fid2 = fopen('../data/gapfill/substrate_usage_rxn.tsv');
format = repmat('%s ',1,4);
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(temp)
    rxnlist(:,i) = temp{i};
end
load('StrainData.mat')

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

% add those rxns for each strain
current_path = pwd;
for i = 2:length(rxnlist(1:end,1))
    strainslist = rxnlist(i,3);
    strainslist = split(strainslist,',');
    strainslist = strrep(strainslist,' ','');
    if strcmp(strainslist,'')
        strainslist = StrianData.strains;
    end
    strainslist=regexprep(strainslist,'_16....$','');
    [~,idx] = ismember(lower(strainslist),lower(StrianData.strains));
    [~,rxnIdx] = ismember(rxnlist(i,1),model_original.rxnMetaNetXID);
    for j = 1:length(idx)
        if idx(j) ~=0
            m = StrianData.strains{idx(j)};
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

