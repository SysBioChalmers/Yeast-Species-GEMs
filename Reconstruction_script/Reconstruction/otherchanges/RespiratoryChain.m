function [newrxns] = RespiratoryChain(model_original,strains,inputpath)
% this function is to plot the repiratory chain complex I and adding back
% the respiratory chain reactions

newrxns = [];
current_path = pwd;
% load the complex I blas result
fid2 = fopen('../../data/physiology/complexI_existence.tsv');
format = '%s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
strain_withoutcomplexI = data{1};
count = cellfun(@str2num,data{2});
idx = find(count < 10);
strain_withoutcomplexI = strain_withoutcomplexI(idx);
[~,rxnMatrix,~] = getprecursorMatrixCobra(model_original,strains,inputpath,'',0);
% add respiratory chain back
rxn = {'r_0438','r_0439','r_0226','r_0773'};
for i = 1:3
    [~,idx] = ismember(rxn(i),model_original.rxns);
    rxnexist = rxnMatrix(:,idx);
    donthave = find(rxnexist==0);
    strains_add = intersect(strains(donthave),strains);
    for j = 1:length(strains_add)
        m = strains_add{j};
        cd(inputpath)
        load([m,'.mat'])
        cd(current_path)
        reducedModel = addrxnBack(reducedModel,model_original,rxn(i),model_original.grRules(idx));
        newrxns = [newrxns;rxn(i),m,{''},{'respiratory'}];
        cd(inputpath)
        save([strains{donthave(j)},'.mat'],'reducedModel')
    end
end
% only gapfill the r_0773 for the species without complexI; THOSE WITH
% COMPLEXI have been fixed at add alternative pathways
for i = 1:length(strain_withoutcomplexI)
    m = strain_withoutcomplexI{i};
    [~,idx] = ismember(lower(m),lower(strains)); % check whether it is in the list of strains
    if idx~=0
         cd(inputpath)
        load([strains{idx},'.mat'])
        cd(current_path)
        [~,idx2] = ismember('r_0773',model_original.rxns);
        reducedModel = addrxnBack(reducedModel,model_original,rxn(4),model_original.grRules(idx2));
        newrxns = [newrxns;rxn(4),m,{''},{'complexI'}];
        cd(inputpath)
        save([m,'.mat'],'reducedModel')
    else
        warning(['no species found for',m])
    end
end

cd(current_path)