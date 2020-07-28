function [model,addedrxn,EnergyResults,MassChargeresults,RedoxResults] = addRxnfromMNX(model,rxnID,comp1,comp2,GPR)
% This function is to add rxns directly from MNX database

EnergyResults     = {};
MassChargeresults = {};
RedoxResults      = {};
addedrxn = [];
load('../../data/databases/MNXmet.mat')
load('../../data/databases/MNXrxn.mat')
if nargin < 5
    GPR = cell(length(rxnID),1);
end

% load MNX reactions
for i = 1:length(rxnID)
    i
    [~,rxnIdx] = ismember(rxnID(i),MNXrxn(:,1)); % index those rxns
    if isequal(rxnIdx,0)
        warning([rxnID{i}, 'may have updated to a new ID in MNX databse, please check and modify to a new ID']);
    else
        %rxn information
        newrxn.ID = rxnID(i);
        newrxn.rxnNames     = MNXrxn(rxnIdx,10);
        newrxn.rxnECNumbers = MNXrxn(rxnIdx,5);
        newrxn.rxnKEGGID    = MNXrxn(rxnIdx,9);
        newrxn.rxnMetaNetXID   = MNXrxn(rxnIdx,1);
        newrxn.direction = MNXrxn(rxnIdx,7);
        % matrix information
        [S, mets]=constructS(MNXrxn(rxnIdx,2));
        matrix.metcoef = S;
        if ~isequal(matrix.metcoef,0)
            temp = split(mets,'@',2);
            matrix.metMetaNetXID = cellstr(temp(:,1));
            comps = cellstr(temp(:,2));
            comps = strrep(comps,'MNXD1',comp1);
            comps = strrep(comps,'MNXD2',comp2);
            matrix.metcompartments = comps;
            
            % find met information
            [~,ID] = ismember(matrix.metMetaNetXID,MNXmet(:,1));
            matrix.metNames =  MNXmet(ID,2);
            matrix.metFormulas      = MNXmet(ID,3);
            matrix.metCharges       = MNXmet(ID,4);
            matrix.metKEGGID        = MNXmet(ID,10);
            matrix.metChEBIID       = MNXmet(ID,11);
            
            % find metname in MNX databse
            matrix.metMetaNetXID(ismember(matrix.metMetaNetXID,'MNXM01')) = {'MNXM1'}; % fix the error for h+ with two IDs
            [~,ID] = ismember(matrix.metMetaNetXID,model.metMetaNetXID);
            metname_temp = split(model.metNames(ID(ID~=0)),' [',2);
            metname_temp = metname_temp(:,1);
            metname_temp = cellstr(metname_temp(:,1));
            matrix.metNames(ID~=0) = metname_temp;
            matrix.metIDs = matrix.metNames;
            
            %change coefficient
            if cell2mat(newrxn.direction) == -1
                matrix.metcoef = matrix.metcoef*-1;
                newrxn.rev = 0;
            elseif cell2mat(newrxn.direction) == 0
                newrxn.rev = 1;
            elseif cell2mat(newrxn.direction) == 1
                newrxn.rev = 0;
            end
            newrxn = rmfield(newrxn,'direction');
            
            %change compartments
            CONValldata = cat(2,model.compNames,model.comps);
            lbracket    = ' [' ;%  space
            llbracket   = '[';
            rbrackets   = ']';
            [m, n]      = size(CONValldata);
            for k = 1:m
                aa = CONValldata(k,1);
                aa = char(aa);
                for j=1:length(matrix.metNames)
                    bb = matrix.metcompartments(j,1);
                    bb = char(bb);
                    if strcmp(bb,aa)
                        matrix.Newcomps(j,1) = CONValldata(k,2);
                    end
                end
            end
            for mm=1:length(matrix.metNames)
                matrix.metnames(mm) = strcat(matrix.metIDs(mm),lbracket,matrix.metcompartments(mm),rbrackets);
                matrix.Newcomps(mm) = strcat(llbracket,matrix.Newcomps(mm),rbrackets);
            end
            
            %mapping mets to model.metnames, get s_ index for new mets
            for j = 1:length(matrix.metnames)
                [~,metindex] = ismember(matrix.metnames(j),model.metNames);
                if metindex ~= 0
                    matrix.mets(j) = model.mets(metindex);
                elseif metindex == 0
                    newID = getNewIndex(model.mets);
                    matrix.mets(j) = strcat('s_',newID,matrix.Newcomps(j));
                    model = addMetabolite(model,char(matrix.mets(j)), ...
                        'metName',matrix.metnames(j));
                end
            end
            
            for f = 1:length(matrix.metnames)
                [~,metID] = ismember(matrix.metNames(f),model.metNames);
                if metID ~= 0
                    model.metFormulas{metID}     = matrix.metFormulas{f};
                    model.metCharges(metID)      = matrix.metCharges(f);
                    model.metKEGGID{metID}       = matrix.metKEGGID{f};
                    model.metChEBIID{metID}      = matrix.metChEBIID{f};
                    model.metMetaNetXID{metID}   = matrix.metMetaNetXID{f};
                end
            end
            
            %add new reactions according to rev ID. Met Coef need to be in the column,
            %not a row. Coef should be double, which was converted at the import
            %section.
            
            newID   = getNewIndex(model.rxns);
            Met = matrix.mets;
            Coef = transpose(matrix.metcoef);
            [model, rxnIDexists] = addReaction(model,...
                ['r_' newID],...
                'reactionName', cell2mat(newrxn.ID),...
                'metaboliteList',Met,...
                'stoichCoeffList',Coef,...
                'reversible',newrxn.rev,...
                'geneRule', GPR{i},...
                'checkDuplicate',1);
            addedrxn = [addedrxn;newrxn.ID];
            [EnergyResults,RedoxResults] = CheckEnergyProduction(model,{['r_' newID]},EnergyResults,RedoxResults);
            [MassChargeresults] = CheckBalanceforSce(model,{['r_' newID]},MassChargeresults);
            
            %add rxn annotation
            [~,idx] = ismember(newrxn.ID,model.rxnNames);
            if idx ~= 0
                model.rxnNames(idx)     = newrxn.rxnNames;
                model.rxnECNumbers(idx) = newrxn.rxnECNumbers;
                model.rxnKEGGID(idx)    =  newrxn.rxnKEGGID;
                model.rxnMetaNetXID(idx)    =  newrxn.rxnMetaNetXID;
                %model.subSystems(rxnID)    =  newrxn.subSystems;
                if isfield(newrxn,'rxnNotes')
                    model.rxnNotes(idx)    =  newrxn.rxnNotes;
                end
            end
            
            clearvars matrix newrxn
        else
            warning(['S is empty, check equation for ',rxnID{i}])
        end
    end
end
end
