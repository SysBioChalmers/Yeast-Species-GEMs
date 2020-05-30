function [newrxns,mets_test] = updateGapRxns(model_original,rxn,met,strains,inputpath,outputpath,type,checkmode)
% This function is to find strains without the rxn but have this rxn in the
% draft model and add that one back and save the model in the filefolder
% specified.

% type is for losse/strict. strict means that the model will only add reactions exist in both draft models from kegg/metacyc. loose will add gap-filling reactions with only one draft models
% checkmode is for whether check the proMatrix or not. Default is 1. If no met is inputted, then it will use all precursors in all pseudo reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = pwd;
if ~exist('checkmode','var') || isempty(checkmode)
    checkmode = 1; % which will test whether the model can produce all precursors; checkmode == 0  means no precursor will be test
end

[proMarix,rxnMatrix,mets_test] = getprecursorMatrixCobra(model_original,strains,inputpath,met,checkmode); % model_original is cobra format panmodel is raven format.
save('allMatrix.mat')
newrxns = [];
if ~isempty(met)
    [~,metIDx] = ismember(met,mets_test);
    if isempty(metIDx)
        warning([met,'does not exist in the biomass, should be met in biomass or use false for no met here'])
    end
end
for j = 1:length(rxn(:,1))
    j
[~,idx] = ismember(rxn(j,1),model_original.rxns);
rxnexist = rxnMatrix(:,idx);
donthave = find(rxnexist==0);
if ~isempty(donthave)
    [result_metacyc,gpr_metacyc] = rxnExistDraft(strains,rxn(j,2),'metacyc');
    [result_kegg,gpr_kegg] = rxnExistDraft(strains,rxn(j,3),'kegg');
    if strcmp(type,'strict')
        result = intersect(find(result_metacyc),find(result_kegg));
    elseif strcmp(type,'loose')
        result = union(find(result_metacyc),find(result_kegg));
    else
        error('should use the term strict or loose for the type')
    end
    straingapfill = intersect(donthave,result);
    if ~isempty(straingapfill)
        if ~isempty(met)
            straingapfill = intersect(straingapfill,find(proMarix(:,metIDx)==0));
        end
        if ~isempty(straingapfill)
        gpr_metacyc_tmp = gpr_metacyc(straingapfill);
        gpr_kegg_tmp = gpr_kegg(straingapfill);
        for i = 1:length(gpr_metacyc_tmp)
            temp_kegg = [];
            temp_metacyc = [];
            gpr_tmp = [];
            if ~isempty(gpr_metacyc_tmp{i})&& ~isempty(gpr_kegg_tmp{i})
                temp_kegg = strsplit(gpr_kegg_tmp{i},' or ');
                temp_metacyc = strsplit(gpr_metacyc_tmp{i},' or ');
                gpr_tmp = unique([temp_kegg,temp_metacyc]);
                gpr(i,1) = {strjoin(gpr_tmp,' or ')};
                Note(i,1) = {'Metacyckegg'};
            elseif ~isempty(gpr_metacyc_tmp{i})&& isempty(gpr_kegg_tmp{i})
                gpr(i,1) = gpr_metacyc_tmp(i);
                Note(i,1) = {'Metacyc'};
            elseif isempty(gpr_metacyc_tmp{i})&& ~isempty(gpr_kegg_tmp{i})
                gpr(i,1) = gpr_kegg_tmp(i);
                Note(i,1) = {'kegg'};
            else
                error('found no gpr')
            end
        end
        newrxns = [newrxns;repmat(rxn(j,1),length(gpr),1),strains(straingapfill),gpr,Note];
        clearvars  Note gpr temp_metacyc gpr_tmp temp_kegg straingapfill donthave rxnexist
        end
    end
end
end
if ~isempty(newrxns)
strains_specific = unique(newrxns(:,2));
for i = 1:length(strains_specific)
    cd(inputpath)
    load([strains_specific{i},'.mat'])
    idx = find(contains(newrxns(:,2),strains_specific(i)));
    for j = 1:length(idx)
     cd(path)
    reducedModel = addrxnBack(reducedModel,model_original,newrxns(idx(j),1),newrxns(idx(j),3));
    end
    cd(outputpath)
    save([strains_specific{i},'.mat'],'reducedModel')
    %saveSSModel(reducedModel,'false');
end
end

cd(path)
end
