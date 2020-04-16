function model = curateKcatValues(model)
%Curate Kcat values
cd GECKO/geckomat/utilities
% proteins = {'O42933';'A0A3F2XYL4';'A0A0X8HUP6';'Q6CWE9';'W0T8T6';'C4QX36';...
%             'A0A1G4MB32';'C5DFP3';'G0V6T0';'A0A0L8RCR8';'Q10283';'I2H574';...
%             'G8BQB6';'C5DS84';'A0A0A8L4H2'};
% value    = (24*60*1e3/1e3*60);
% Kcats    = [70.9;value;value;value;70.29;78.3;70.29;value;value;value;...
%             value;value;value;value;value];
% for i=1:length(proteins)
%     prot = proteins{i};
%     Kcat = Kcats(i);
%     model = assignKcat(model,prot,Kcat);
% end
% %Curate Kcats by rxnNames too
 rxns = {'hydroxymethylglutaryl CoA synthase' ...
         'hydroxymethylglutaryl CoA reductase'};%'tryptophan synthase (indoleglycerol phosphate)' ...
%         'IMP dehydrogenase'};
kcats = [33.3333 ... %blattella germanica
         24]; %Haloferax volcanii
%         2.46]; %imp	eremothecium gossypii//*//*
for i=1:length(rxns)  
    rxn  = rxns{i};
    kcat = kcats(i);
    indexes = find(contains(model.rxnNames,rxn));
    for j=1:length(indexes)
        index = indexes(j);
        mets  = find(model.S(:,index));
        prots = mets(find(contains(model.mets(mets),'prot_')));
        if ~isempty(prots) & ~isempty(index)
            model.S(prots,index) = -1/(kcat*3600);
        end
    end
end
cd ../../..
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = assignKcat(model,protein,Kcat)
[~,rxnIdx,~,~] = getKcat(model,protein);
prot = find(contains(model.metNames,protein));
newValue = Kcat; %1/s
model.S(prot,rxnIdx) = -1/(newValue*3600);
end