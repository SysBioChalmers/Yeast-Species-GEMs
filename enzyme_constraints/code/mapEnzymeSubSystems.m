function subSystems = mapEnzymeSubSystems(enzymes,model)
%mapEnzymeSubSystems
%
% Function that gets a list of enzymes and serach all the subSystems/pathways
% in which it is involved according to the provided model structure.
%
%   model      An ecModel structure.
%   enzymes    (cell) List of enzyme IDs (uniprot) to search in the model.
%
%   subSystems (cell) Array of strings with all the subsytems for each of
%              the queried enzymes.
%
%   usage: subSystems = mapEnzymeSubSystems(enzymes,model)
%
% Last modified.  Ivan Domenzain 2019-04-16

subSystems = {};
for i=1:length(enzymes)
    enzyme_i = ['prot_' enzymes{i}];
    metPos   = find(strcmpi(model.metNames,enzyme_i));
    rxnsProt = find(model.S(metPos,:));
    rxnsProt = rxnsProt(1:end-1);
    subSystem = {};
    for j=1:length(rxnsProt)
        index = rxnsProt(j);
        if isempty(model.subSystems{index})
            str = ' ';
        else
            str = strjoin(model.subSystems{index},' // ');
        end
        subSystem = [subSystem; {str}];
    end
    subSystem  = strjoin(subSystem,' // ');
    subSystems = [subSystems; {subSystem}];
end
end