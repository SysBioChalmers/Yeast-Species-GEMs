function [model,modifications] = manualModifications(model)
%
% model = manualModifications(model)
% 
% Ivan Domenzain.      Last edited: 2020-03-09

modifications{1} = [];
modifications{2} = [];

% Remove repeated reactions (2017-01-16):
rem_rxn = false(size(model.rxns));
for i = 1:length(model.rxns)-1
    for j = i+1:length(model.rxns)
        if isequal(model.S(:,i),model.S(:,j)) && model.lb(i) == model.lb(j) && ...
                model.ub(i) == model.ub(j)
            rem_rxn(j) = true;
            disp(['Removing repeated rxn: ' model.rxns{i} ' & ' model.rxns{j}])
        end
    end
end
model = removeRxns(model,model.rxns(rem_rxn));
% Merge arm reactions to reactions with only one isozyme (2017-01-17):
arm_pos = zeros(size(model.rxns));
p       = 0;
for i = 1:length(model.rxns)
    rxn_id = model.rxns{i};
    if contains(rxn_id,'arm_')
        rxn_code  = rxn_id(5:end);
        k         = 0;
        for j = 1:length(model.rxns)
            if contains(model.rxns{j},[rxn_code 'No'])
                k   = k + 1;
                pos = j;
            end
        end
        if k == 1
            %Condense both reactions in one:
            new_id     = model.rxns{pos};
            new_name   = model.rxnNames{pos};
            stoich     = model.S(:,i) + model.S(:,pos);
            model      = addReaction(model,{new_id,new_name},model.mets,stoich,true,0,1000);
            p          = p + 1;
            arm_pos(p) = i;
            disp(['Merging reactions: ' model.rxns{i} ' & ' model.rxns{pos}])
        end
    end
end
% Remove saved arm reactions:
model = removeRxns(model,model.rxns(arm_pos(1:p)));
%Change gene rules:
if isfield(model,'rules')
    for i = 1:length(model.rules)
        if ~isempty(model.rules{i})
            %Change gene ids:
            model.rules{i} = strrep(model.rules{i},'x(','');
            model.rules{i} = strrep(model.rules{i},')','');
            model.rules{i} = model.genes{str2double(model.rules{i})};
        end
    end
end
% Remove unused enzymes after manual curation (2017-01-16):
rem_enz = false(size(model.enzymes));
for i = 1:length(model.enzymes)
    pos_met = strcmp(model.mets,['prot_' model.enzymes{i}]);
    if sum(model.S(pos_met,:)~=0) == 1
        rem_enz(i) = true;
    end
end
rem_enz = model.enzymes(rem_enz);
for i = 1:length(rem_enz)
    model = deleteProtein(model,rem_enz{i});
    disp(['Removing unused protein: ' rem_enz{i}])
end
% Block O2 and glucose production (for avoiding multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;
% Map the index of the modified Kcat values to the new model (after rxns
% removals).
modifications = mapModifiedRxns(modifications,model);
%Curate Kcat values
cd ../utilities
proteins = {'O42933' 'P46971' 'P36143' 'P52867' 'C4R045' 'A0A3F2XYL4' 'C5DS84' 'W0TFZ8'};
Kcats    = [20 20 13.92 20 20 20.1 20.1 29.9];
for i=1:length(proteins)
    prot = proteins{i};
    Kcat = Kcats(i);
    model = assignKcat(model,prot,Kcat);
end
cd ../change_model
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = assignKcat(model,protein,Kcat)
[~,rxnIdx,~,~] = getKcat(model,protein);
prot = find(contains(model.metNames,protein));
newValue = Kcat; %1/s
model.S(prot,rxnIdx) = -1/(newValue*3600);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modified = mapModifiedRxns(modifications,model)
modified = [];
for i=1:length(modifications{1})
    rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
    str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
    modified = [modified; str];
end
end