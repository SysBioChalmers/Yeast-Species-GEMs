function model = flexibilize_bottleNecks(model)
cd GECKO/geckomat/utilities
proteins = {'O42933' 'P46971' 'P36143' 'P52867' 'C4R045' 'C5DS84'};%  'A0A3F2XYL4'  'W0TFZ8'};
Kcats    = [20 20 13.92 20 20 3160];% 0.617 0.617 0.617];
common = intersect(model.enzymes,proteins);
if ~isempty(common)
    for i=1:length(common)
        index = strcmpi(proteins,common(i));
        kcat = Kcats(index);
        [~,rxnIdx] = getKcat(model,common{i});
        protIdx = find(strcmpi(model.metNames,['prot_' common{i}]));
        model.S(protIdx,rxnIdx) = -1/(3600*kcat);
    end
end
cd ../../..
end   
