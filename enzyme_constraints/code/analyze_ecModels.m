%analyze_ecModels
%
% Ivan Domenzain. 2020-04-15
%
current = pwd;
%Clone GECKO repository
%git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
%git('pull')
%Locate the correct branch
%git('stash')
%git('checkout feat-add_utilities') 
clc
cd ..
%load organism and model specific parameters
fID         = fopen('../data/yeasts_parameters.txt');
yeastsParam = textscan(fID,'%s %s %s %s %f %f %f','Delimiter','\t','HeaderLines',1);
%load phenotypes information
fID      = fopen('../../Reconstruction_script/data/physiology/343_phenotype_clade.txt');
phenData = textscan(fID,'%s %s %s','Delimiter','\t','HeaderLines',1);
phenData = table(phenData{1},phenData{3},'VariableNames',{'name' 'phenotype'});
phenData.name      = strrep(phenData.name,'"','');
phenData.phenotype = strrep(phenData.phenotype,'"','');
biomass_phen       = readtable('../../Reconstruction_script/data/physiology/biomass_type.txt');
%Retrieve ecModel names
fileNames = dir('../ecModels');
originalModels = dir('../models');
modelsData = table();
for i=1:length(originalModels)
    cd(current)
    file = originalModels(i).name;
    if contains(file,'.mat')
        modelFile = file(1:end-4);
        ecModelName = ['ec_' modelFile '_GEM'];
        disp(modelFile)
        %get phenotype
        nameStr = strrep(modelFile,'_',' ');
        phenotype = phenData.phenotype{strcmpi(phenData.name,nameStr)};
        phenotype = strsplit(phenotype,',');
        %Load original model
        load(['../models/' file])
        %Load ecModel
        load(['../ecModels/' ecModelName '/ecModel_batch.mat'])
        %Transfer org specific parameters to GECKO
        transferParameters(yeastsParam,reducedModel,modelFile);      
        cd ..
        %Get classification od matched enzymes
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
        forEx     = find(strcmpi(ecModel_batch.rxnNames,'formate exchange'));
        mEtExc    = find(strcmpi(ecModel_batch.rxnNames,'methanol exchange (reversible)'));
        Ppool    = find(strcmpi(ecModel_batch.rxnNames,'prot_pool_exchange'));
        %Get a FBA solution
        solution  = solveLP(ecModel_batch);
        met_yield = []; 
        if ~isempty(solution.x)
            printFluxes(ecModel_batch, solution.x, true, 1E-5)
            gRate    = solution.x(obj);
            load('GECKO/geckomat/parameters.mat')
            %Introduce curation of problematic Kcats
            ecModel_batch = curateKcatValues(ecModel_batch);
            %Fit sigma parameter and adjust protein pool upper bound
            ecModel_batch = rescaleBioMassComposition(ecModel_batch,phenotype{1},biomass_phen);
            cd ../kcat_sensitivity_analysis
            %Get U_max
            solution = solveLP(ecModel_batch);
            OptSigma = sigmaFitter(ecModel_batch,0.46,gRate,0.5);
            ecModel_batch.ub(end) = 0.46*OptSigma*0.5;
            %Save calibrated ecModels
            save(['../../../../ecModels/ec_' modelFile '_GEM/ecModel_batch_curated.mat'],'ecModel_batch')
            %Constrain acetate production
            %ecModel_batch.ub(acEx) = 0.6;
            %ecModel_batch.ub(forEx) = 0;
            %Get exchange fluxes from a new solution
            solution = solveLP(ecModel_batch,1);
            GUR      = solution.x(cSource);
            acExc    = solution.x(acEx);
            %EtExc    = solution.x(EtOH);
            %Get the top used enzymes
            T = topUsedEnzymes(solution.x,ecModel_batch,{'glc'},ecModelName,false);
            [lim_growth,breakFlag] = findTopLimitations(ecModel_batch,[]);
            %check ethanol production capabilities
            temp     = setParam(ecModel_batch,'obj',Ppool,-1);
            temp     = setParam(temp,'lb',obj,0.99*solution.x(obj));
            temp     = setParam(temp,'ub',obj,0.6);
            %temp     = setParam(temp,'ub',obj,0);
            solution = solveLP(temp,1);
            EtExc    = solution.x(EtOH);            
            temp     = setParam(temp,'lb',Ppool,0.99999999*solution.x(Ppool));            
            temp     = setParam(temp,'obj',EtOH,1);
            [lim_EtOH,~] = findTopLimitations(temp,[]);
            %Check growth on methanol
            if ~isempty(mEtExc)
                temp = setParam(ecModel_batch,'obj',obj,1);
                temp = setParam(temp,'ub',cSource,0);
                temp = setParam(temp,'ub',mEtExc,1000);
                solution  = solveLP(temp,1);
                if ~isempty(solution.x)
                    met_yield = solution.x(obj)/(solution.x(mEtExc)*0.03204);
                end
            end
        else
            gRate = 0;
            GUR   = 0;
            EtExc  = 0;
        end
        cd ../utilities
        if isempty(met_yield)
            met_yield = 0;
        end
        f_factor = OptSigma;
        [kcat,rxnIdx,rxnName,MW] = getKcat(ecModel_batch,T.prots_glc{1});
        modelsData = [modelsData; {modelFile phenotype{1} nRxns nMets nGnes nRxns_ec ...
                      nMets_ec nEnz_ec isoEnzymes Promiscuous Complexes ...
                      RxnWithKcat gRate GUR acExc EtExc met_yield yeastParam.GUR yeastParam.EtOH ...
                      T.prots_glc{1} T.Usages_glc(1) kcat(1) rxnName(1) MW(1) ...
                      lim_growth{1} lim_growth{4} lim_growth{5} lim_growth{6} ...
                      lim_EtOH{1} lim_EtOH{5} lim_EtOH{6} f_factor}];
        close all
        clc
        save(['../../../../ecModels/ec_' modelFile '_GEM/ecModel_batch.mat'],'ecModel_batch')
    end
end
modelsData.Properties.VariableNames ={'model' 'phenotype' 'nRxns' 'nMets' 'nGenes' 'ec_nRxns' 'ec_nMets' 'ec_nEnz' ...
            'isoenzymes' 'promiscuous' 'complexes' 'RxnWithKcat' 'gRate' ...
            'GUR' 'acExc' 'EtExc' 'meth_yield' 'GUR_exp' 'EtOH_exp' 'topProt' 'topUsage' 'kcat' ...
            'rxnName' 'MW' 'limEnz' 'limKcat' 'gCC' 'limRxn' 'limEnz_EtOH' 'EtOH_CC' 'limRxn_EtOH' 'f_factor'};
%Save results table
mkdir('../../../../results')
writetable(modelsData,'../../../../results/ecModels_metrics.txt','Delimiter','\t','QuoteStrings',false)
        