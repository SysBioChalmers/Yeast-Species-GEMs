current = pwd;
%Retrieve ecModel names
fileNames = dir('../ecModels');
originalModels = dir('../models');

for i=1:length(originalModels)
    cd(current)
    file = originalModels(i).name;
    if contains(file,'.mat')
        modelFile = file(1:end-4);
        load(['../models/' file])
        %load(['../ecModels/ec_' modelFile '_GEM/ecModel_batch.mat'])
        t = readtable(['../results/gCC/' modelFile '_limGrowth.txt']);
        [member,iA] = ismember(t.enzGenes,reducedModel.genes);
        t.enz_panIDs = reducedModel.proteins(iA);
        writetable(t,['../results/gCC/' modelFile '_limGrowth.txt'],'delimiter','\t','QuoteStrings',false)
    end
end