function curate_ecModels(initModel,finalModel,resultsFile,fittingOption,getECC,fixEtOH)
%curate_ecModels
%
% Last modified. Ivan Domenzain. 2020-05-27
%
current = pwd;
%Clone GECKO repository
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
%clc
cd ..
%load organism and model specific parameters
fID         = fopen('../data/yeasts_parameters.txt');
yeastsParam = textscan(fID,'%s %s %s %s %f %f %f','Delimiter','\t','HeaderLines',1);
%load phenotypes information
fID      = fopen('../../Reconstruction_script/data/physiology/343_phenotype_clade.tsv');
phenData = textscan(fID,'%s %s %s','Delimiter','\t','HeaderLines',1);
phenData = table(phenData{1},phenData{3},'VariableNames',{'name' 'phenotype'});
phenData.name      = strrep(phenData.name,'"','');
phenData.name      = strrep(phenData.name,'_',' ');
phenData.phenotype = strrep(phenData.phenotype,'"','');
biomass_phen       = readtable('../../Reconstruction_script/data/physiology/biomass_type.txt');
%Retrieve ecModel names
fileNames = dir('../ecModels');
originalModels = dir('../models');
modelsData = table();
mkdir('../results/fluxDist')
mkdir('../results/gCC')
cd(current)
for i=1:length(originalModels)
    file = originalModels(i).name;
    if contains(file,'.mat')
        modelFile = file(1:end-4);
        ecModelName = ['ec_' modelFile '_GEM'];
        disp(['***** ' modelFile])
        %get phenotype
        nameStr = strrep(modelFile,'_',' ');
        position = find(strcmpi(phenData.name,nameStr));
        phenotype = phenData.phenotype{position};
        phenotype = strsplit(phenotype,',');
        %Load original model
        load(['../models/' file])
        %Load ecModel
        load(['../ecModels/' ecModelName '/' initModel '.mat'])
        %Correct Kcats for identified bottlenecks
        ecModel_batch = flexibilize_bottleNecks(ecModel_batch);
        %Transfer org specific parameters to GECKO
        transferParameters(yeastsParam,reducedModel,modelFile);
        cd ..
        %Get classification of matched enzymes
        [isoEnzymes,Promiscuous,Complexes,RxnWithKcat] = rxnCounter(ecModel_batch,reducedModel);
        nRxns    = length(reducedModel.rxns);
        nMets    = length(reducedModel.mets);
        nGnes    = length(reducedModel.genes);
        nRxns_ec = length(ecModel_batch.rxns);
        nMets_ec = length(ecModel_batch.mets)-length(ecModel_batch.enzymes);
        nEnz_ec  = length(ecModel_batch.enzymes);
        %Get relevant indexes
        obj      = find(contains(ecModel_batch.rxnNames,'biomass pseudoreaction'));
        cSource  = find(contains(ecModel_batch.rxnNames,'D-glucose exchange (reversible)'));
        EtOH     = find(strcmpi(ecModel_batch.rxnNames,'ethanol exchange'));
        acEx     = find(strcmpi(ecModel_batch.rxnNames,'acetate exchange'));
        mEtEx    = find(strcmpi(ecModel_batch.rxnNames,'methanol exchange (reversible)'));
        Ppool    = find(strcmpi(ecModel_batch.rxnNames,'prot_pool_exchange'));
        MetExc   = 0;
        MetGrate = 0;
        %Get exp miu_max
        load('GECKO/geckomat/parameters.mat')
        gRate = yeastParam.gR_exp;
        disp(['The experimental growth rate is: ' num2str(gRate)])
        %rescale biomass components
        ecModel_batch = rescaleBioMassComposition(ecModel_batch,phenotype{1},biomass_phen);
        ecModel_batch = setParam(ecModel_batch,'obj',obj,1);
        %Find model phenotype info
        index = strcmpi(biomass_phen.Properties.VariableNames,phenotype);
        Ptot  = table2array(biomass_phen(1,index));
        %Fit sigma parameter in order to minimize GUR prediction error
        cd ../../../specific_scripts
        ecModel_batch.lb(obj) = 0.95*gRate;
        ecModel_batch.ub(obj) = 1.05*gRate;
        switch fittingOption 
            case 1
                [Opt_f,error] = sigmaFitter(ecModel_batch,Ptot,yeastParam.GUR,0.5,cSource);
                ecModel_batch.ub(end) = Ptot*Opt_f*0.5;
                error = error/100;
                ecModel_batch.lb(cSource) = (1-1.05*error)*yeastParam.GUR;
                ecModel_batch.ub(cSource) = (1+1.05*error)*yeastParam.GUR;
            case 2
                [Opt_f,error] = sigmaFitter(ecModel_batch,Ptot,yeastParam.gR_exp,0.5,obj);   
                ecModel_batch.ub(end) = Ptot*Opt_f*0.5;
                error = error/100;
                ecModel_batch.lb(obj) = (1-1.05*error)*yeastParam.gR_exp;
                ecModel_batch.ub(obj) = (1+1.05*error)*yeastParam.gR_exp;
            case 3
                [Opt_f,error] = sigmaFitter(ecModel_batch,Ptot,yeastParam.EtOH,0.5,EtOH);   
                ecModel_batch.ub(end) = Ptot*Opt_f*0.5;
                error = error/100;
                ecModel_batch.lb(EtOH) = (1-1.05*error)*yeastParam.EtOH;
                ecModel_batch.ub(EtOH) = (1+1.05*error)*yeastParam.EtOH;

            otherwise
                Opt_f = 0;
        end
        %Get exchange fluxes
        solution = solveLP(ecModel_batch,1);
        GUR      = solution.x(cSource);
        acExc    = solution.x(acEx);
        yeastParam.EtOH;
        EtExc    = solution.x(EtOH);
        %Calculate error in biomass yield  
        gSim =solution.x(obj);
        bioYield = gSim/(GUR*0.18);
        yieldExp = gRate/(yeastParam.GUR*0.18);
        bioError = (bioYield-yieldExp)/yieldExp;
        %disp(['The biomass production prediction is: ' num2str(solution.x(obj))])
        disp(['The error in the biomass yield is: ' num2str(bioError*100) '%'])
        cd ..
        %Get rxns table
        varNames   = {'rxns' 'rxnNames' 'formulas' 'flux' 'grRules' 'subSystems'};
        rxnIndxs   = find(~contains(ecModel_batch.rxnNames,'draw_prot_'));
        enzIndexes = find(contains(ecModel_batch.rxnNames,'draw_prot_'));
        rxnsTable  = getSubsetTable(rxnIndxs,ecModel_batch,solution.x,varNames,false);
        newFile    = ['../results/fluxDist/ec_' modelFile '_rxnsTable.txt'];
        writetable(rxnsTable,newFile,'Delimiter', '\t','QuoteStrings',false);
        %Get absolute enzyme usages table
        varNames     = {'enzymes' 'shortNames' 'abs_usage' 'grRules' 'subSystems'};
        enzTable_abs = getSubsetTable(enzIndexes,ecModel_batch,solution.x,varNames,true);
        newFile      = ['../results/fluxDist/ec_' modelFile '_enzTable.txt'];
        writetable(enzTable_abs,newFile,'Delimiter', '\t','QuoteStrings',false);
        cd GECKO/geckomat/kcat_sensitivity_analysis
        %Get the top used enzymes and perform Kcat sensitivity analysis on
        %gRate
        T = topUsedEnzymes(solution.x,ecModel_batch,{'glc'},ecModelName,false);
        if getECC
            unconst_model = ecModel_batch;
            unconst_model.ub(obj) = 1000;
            unconst_model.lb(obj) = 0;
            unconst_model.lb(cSource) = 0;
            unconst_model.ub(cSource) = 1000;
            unconst_model.ub(EtOH) = 1000;
            unconst_model.lb(EtOH) = 0;
            [lim_growth,~] = findTopLimitationsAll(unconst_model,[],obj,1.001);
            if ~isempty(lim_growth)
                %temp           = table(lim_growth{1},lim_growth{2},lim_growth{3},lim_growth{4},lim_growth{5},lim_growth{6});
                writetable(lim_growth,['../../../../results/gCC/' modelFile '_limGrowth.txt'],'Delimiter','\t','QuoteStrings',false)
            else 
                disp('No limitations were found')
            end
            ecModel_batch.ub(obj) = 1.01*solution.x(obj);
        end
        %check ethanol production capabilities
        tempM    = setParam(ecModel_batch,'obj',EtOH,1);
        tempM    = setParam(tempM,'lb',obj,0.9999*solution.x(obj));
        solution = solveLP(tempM,1);
        EtMax    = solution.x(EtOH);
        if fixEtOH
            if EtMax>0
                tempM = ecModel_batch;
                if EtMax>=yeastParam.EtOH
                    tempM.lb(EtOH) = 0.95*yeastParam.EtOH;
                else
                    tempM.lb(EtOH) = 0.95*EtMax;
                end
                etProd = max(EtExc,yeastParam.EtOH);
                tempM.ub(EtOH) = 1.01*etProd;
            end
            %Obtain definite exchange fluxes
            solution = solveLP(tempM,1);
            if ~isempty(solution.x)
                ecModel_batch = tempM;
                GUR      = solution.x(cSource);
                acExc    = solution.x(acEx);
                EtExc    = solution.x(EtOH);
            end
            disp(['The experimental EtOH production is: ' num2str(yeastParam.EtOH)])
            disp(['The predicted EtOH production is: ' num2str(EtExc)])
            disp(['The predicted max EtOH production is: ' num2str(EtMax)])
            %Check growth on methanol
            if ~isempty(mEtEx)
                tempM = setParam(ecModel_batch,'obj',obj,1);
                tempM = setParam(tempM,'ub',cSource,0);
                tempM = setParam(tempM,'ub',mEtEx,1000);
                solution  = solveLP(tempM,1);
                if ~isempty(solution.x)
                    MetExc = solution.x(mEtEx);%solution.x(obj)/(solution.x(mEtExc)*0.03204);
                    MetGrate = solution.x(obj);
                end
            end
        end
        cd ../utilities    
        %Get rxnName for the top limiting Kcat
        [kcat,~,rxnName,MW] = getKcat(ecModel_batch,T.prots_glc{1});      
        modelsData = [modelsData; {modelFile phenotype{1} nRxns nMets nGnes nRxns_ec ...
                      nMets_ec nEnz_ec isoEnzymes Promiscuous Complexes RxnWithKcat ...
                      gSim GUR acExc EtExc EtMax gRate yeastParam.GUR ...
                      yeastParam.EtOH T.prots_glc{1} T.Usages_glc(1) kcat(1) rxnName(1) MW(1) Opt_f}];
         close all
         disp(' ')
         cd(current)
         %Save calibrated ecModels
         save(['../ecModels/ec_' modelFile '_GEM/' finalModel '.mat'],'ecModel_batch')  
         disp(num2str(yeastParam.EtOH))
         disp(num2str(ecModel_batch.lb(EtOH)))
    end
end
modelsData.Properties.VariableNames = {'model' 'phenotype' 'nRxns' 'nMets' 'nGenes' 'ec_nRxns' ...
                                      'ec_nMets' 'ec_nEnz' 'isoenzymes' 'promiscuous' 'complexes' 'RxnWithKcat' ...
                                      'gRate' 'GUR' 'acExc' 'EtExc' 'EtMax' 'gRate_exp' 'GUR_exp' ...
                                      'EtOH_exp' 'topProt' 'topUsage' 'kcat' 'rxnName' 'MW' 'f_factor'};
%Save results table
cd (current)
mkdir('../results')
writetable(modelsData,['../results/' resultsFile '.txt'],'Delimiter','\t','QuoteStrings',false)
end