function model = curateKcatValues(model)
%Curate Kcat values
cd GECKO/geckomat/utilities
%Curate Kcats by rxnNames
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