current = pwd;
%Rescale biomass composition, fit sigma for growth
curate_ecModels('ecModel_batch','ecModel_modified','ecModels_metrics',2,false,true)
%%Curate Ethanol consumption
clc
cd(current)
%Open model metrics file
metrics    = readtable('../results/ecModels_metrics.txt');
modelNames = metrics.model;
EtOH_exc  = metrics.EtExc;
EtOH_exc(EtOH_exc==0) = 1E-8;
EtOH_exp = metrics.EtOH_exp;
EtOH_exp(EtOH_exp==0) = 1E-8;1
Fchange = EtOH_exc./EtOH_exp;
misspredicted = zeros(length(Fchange),1);
misspredicted(find(Fchange<0.5)) = 1;
misspredicted(find(Fchange>1.5)) = 2;
%Curate ethanol production rate (if needed)
for i=1:length(modelNames)
    modelName = ['ec_' modelNames{i} '_GEM'];
    disp(['******* ' modelName ' *******'])
    disp(['Experimental Ethanol production: ' num2str(EtOH_exp(i))])
    disp(['Predicted Ethanol production: ' num2str(EtOH_exc(i))])
    load(['../ecModels/' modelName '/ecModel_modified.mat'])
    gIndx    = find(strcmpi(ecModel_batch.rxnNames,'biomass pseudoreaction'));
    glcIndx  = find(strcmpi(ecModel_batch.rxnNames,'D-glucose exchange (reversible)'));
    EtOHIndx = find(strcmpi(ecModel_batch.rxnNames,'ethanol exchange'));
    ecModel_batch = setParam(ecModel_batch,'obj',gIndx,1);
    newEtOH = metrics.EtExc;
    if misspredicted(i)>0
        exp_gRate  = metrics.gRate(i);
        %load(['../ecModels/' modelName '/ecModel_modified.mat'])
        %disp(['Predicted Ethanol prod rate: ' num2str(ecModel_batch.ub(EtOHIndx))])
        ecModel_batch.ub(glcIndx) = 1000;
        if misspredicted(i) ==1
            %Underpredictions
            ecModel_batch.ub(EtOHIndx) = 1000;
            [ecModel_batch,newEtOH] = iterative_curation(ecModel_batch,EtOH_exp(i),true,EtOHIndx,gIndx,EtOHIndx,10,20);
        else
            %Overpredictions
        	[ecModel_batch,newEtOH] = iterative_curation(ecModel_batch,EtOH_exp(i),false,EtOHIndx,gIndx,EtOHIndx,0.1,20);
        end
        cd(current)
        ecModel_batch = setParam(ecModel_batch,'obj',gIndx,1);
    end
    
     %underpredicted EtOH exchange
        LB = 0.75;
        UB = 1.25;
        if newEtOH<=LB*EtOH_exp(i)
            ecModel_batch.lb(EtOHIndx) = 0.9999*newEtOH;
            ecModel_batch.ub(EtOHIndx) = 1.05*EtOH_exp(i);
        elseif newEtOH>=UB*EtOH_exp(i)
            ecModel_batch.ub(EtOHIndx) = 1.05*EtOH_exp(i);
            ecModel_batch.lb(EtOHIndx) = 0.95*EtOH_exp(i);
        else
            ecModel_batch.ub(EtOHIndx) = UB*EtOH_exp(i);
            ecModel_batch.lb(EtOHIndx) = LB*EtOH_exp(i);
        end
    save(['../ecModels/' modelName '/ecModel_batch_curated_EtOH.mat'],'ecModel_batch')
    disp(' ')
end
%Readjust 
cd(current)
curate_ecModels('ecModel_batch_curated_EtOH','ecModel_batch_curated_EtOH','ecModels_metrics_curated_EtOH',2,false,true)

%%Curate glucose consumption
cd(current)
%Open model metrics file
metrics    = readtable('../results/ecModels_metrics_curated_EtOH.txt');
%metrics    = metrics([8,12],:);
modelNames = metrics.model;
GUR_sim = metrics.GUR;
GUR_exp = metrics.GUR_exp;
Fchange = GUR_sim./GUR_exp;
misspredicted = zeros(length(Fchange),1);
misspredicted(find(Fchange<0.75)) = 1;
misspredicted(find(Fchange>1.25)) = 2;
for i=1:length(modelNames)
    modelName = ['ec_' modelNames{i} '_GEM'];
    disp(['******* ' modelName ' *******'])
    disp(['Experimental GUR: ' num2str(GUR_exp(i))])
    disp(['Predicted GUR: ' num2str(GUR_sim(i))])
    load(['../ecModels/' modelName '/ecModel_batch_curated_EtOH.mat'])
    gIndx    = find(strcmpi(ecModel_batch.rxnNames,'biomass pseudoreaction'));
    glcIndx  = find(strcmpi(ecModel_batch.rxnNames,'D-glucose exchange (reversible)'));
    EtOHIndx = find(strcmpi(ecModel_batch.rxnNames,'ethanol exchange'));
    if misspredicted(i)>0
        exp_gRate  = metrics.gRate(i);       
        disp(['Predicted Ethanol prod rate: ' num2str(ecModel_batch.ub(EtOHIndx))])
        if misspredicted(i) ==1
            %Underpredictions
            [ecModel_batch,newGUR] = iterative_curation(ecModel_batch,GUR_exp(i),true,glcIndx,gIndx,gIndx,0.1,20);
        else
            %Overpredictions
        	[ecModel_batch,newGUR] = iterative_curation(ecModel_batch,GUR_exp(i),false,glcIndx,gIndx,gIndx,10,20);
        end      
        %underpredicted GUR
        LB = 0.85;
        UB = 1.15;
        if newGUR<=LB*GUR_exp(i)
            ecModel_batch.lb(glcIndx) = 0.9999*newGUR;
            ecModel_batch.ub(glcIndx) = 1.05*GUR_exp(i);
        elseif newGUR>=UB*GUR_exp(i)
            ecModel_batch.ub(glcIndx) = 1.0001*newGUR;
            ecModel_batch.lb(glcIndx) = 0.95*GUR_exp(i);
        else
            ecModel_batch.ub(glcIndx) = UB*GUR_exp(i);
            ecModel_batch.lb(glcIndx) = LB*GUR_exp(i);
        end
        cd(current)
        save(['../ecModels/' modelName '/ecModel_batch_curated.mat'],'ecModel_batch')
    end
    disp(' ')
end
%Fit sigma for fitted GUR, 
cd(current)
curate_ecModels('ecModel_batch_curated','ecModel_batch_improved','ecModels_metrics_improved',1,true,true)
% 


