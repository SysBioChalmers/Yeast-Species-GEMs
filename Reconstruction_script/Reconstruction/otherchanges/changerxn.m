function model = changerxn(model,rxnID,rxnformula)


[~,idx] = ismember(rxnID,model.rxns); 
% This function is to change new rxns
rxnformula = strrep(rxnformula,' [','[');
[metaboliteList, stoichCoeffList, revFlag] = parseRxnFormula(rxnformula);
metaboliteList = strrep(metaboliteList,'[',' [');
metaboliteList = strrep(metaboliteList,'&',' ');
comps = split(metaboliteList', ' [');
comps = comps(:,2);
comps = strrep(comps,']','');
CONValldata = cat(2,model.compNames,model.comps);
[~,b] = ismember(comps,CONValldata(:,1));
comps = CONValldata(b,2);

%mapping mets to model.metnames, get s_ index for new mets
for j = 1:length(metaboliteList)
    [~,metindex] = ismember(metaboliteList(j),model.metNames);
    if metindex ~= 0
        mets(j) = model.mets(metindex);
    elseif metindex == 0
        newID = getNewIndex(model.mets);
        mets(j) = strcat('s_',newID,'[',comps(j),']');
        model = addMetabolite(model,char(mets(j)), ...
                              'metName',metaboliteList{j});
    end
end

    
[model, rxnIDexists] = addReaction(model,...
                                    rxnID,...
                                    'reactionName', model.rxnNames{idx},...
                                    'metaboliteList',mets,...
                                    'stoichCoeffList',stoichCoeffList,...
                                    'reversible',revFlag,...
                                    'checkDuplicate',1);
 
end
