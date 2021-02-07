%Curate glucose consumption
current = pwd;
%Open model metrics file
metrics    = readtable('../results/ecModels_metrics.txt');
modelNames = metrics.model;
bioYield = metrics.gRate./(0.18*metrics.GUR);
bioY_exp = metrics.gRate./(0.18*metrics.GUR_exp);
Fchange = bioYield./bioY_exp;
misspredicted = zeros(length(Fchange),1);
misspredicted(find(Fchange<0.5)) = 1;
misspredicted(find(Fchange>1.5)) = 2;

for i=1:length(modelNames)
    modelName = ['ec_' modelNames{i} '_GEM'];
    disp(['******* ' modelName ' *******'])
    disp(['Experimental Biomass yield: ' num2str(bioY_exp(i))])
    disp(['Biomass yield prediction: ' num2str(bioYield(i))])
    if misspredicted(i)>0
        exp_gRate  = metrics.gRate(i);
        load(['../ecModels/' modelName '/ecModel_batch_modified.mat'])
        gIndx   = find(strcmpi(ecModel_batch.rxnNames,'biomass pseudoreaction'));
        glcIndx = find(strcmpi(ecModel_batch.rxnNames,'D-glucose exchange (reversible)'));
        if misspredicted(i) ==1
            %Underpredictions
            ecModel_batch = iterative_curation(ecModel_batch,exp_gRate,bioY_exp(i),false,glcIndx,gIndx,10,0.18);
        else
            %Overpredictions
        	ecModel_batch = iterative_curation(ecModel_batch,exp_gRate,bioY_exp(i),true,glcIndx,gIndx,10,0.18);
        end
        ecModel_batch = setParam(ecModel_batch,'lb',gIndx,0);
        ecModel_batch = setParam(ecModel_batch,'obj',gIndx,1);
        cd(current)
        save(['../ecModels/' modelName '/ecModel_batch_curated.mat'],'ecModel_batch')
    end
    disp(' ')
end
