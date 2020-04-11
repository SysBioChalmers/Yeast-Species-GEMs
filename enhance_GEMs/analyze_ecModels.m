%analyze_ecModels
%
% Ivan Domenzain. 2020-03-11
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
fID         = fopen('../ComplementaryData/yeasts_parameters.txt');
yeastsParam = textscan(fID,'%s %s %s %s %f %f %f','Delimiter','\t','HeaderLines',1);
%load phenotypes information
fID      = fopen('../Reconstruction_script/data/physiology/343_phenotype_clade.txt');
phenData = textscan(fID,'%s %s %s','Delimiter','\t','HeaderLines',1);
phenData = table(phenData{1},phenData{3},'VariableNames',{'name' 'phenotype'});
phenData.name = strrep(phenData.name,'"','');
phenData.phenotype = strrep(phenData.phenotype,'"','');
biomass_phen = readtable('../Reconstruction_script/data/physiology/biomass_type.txt');
%Retrieve ecModel names
fileNames = dir('../ecModels');
originalModels = dir('models');
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
        load(['models/' file])
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
        %Get a FBA solution
        solution = solveLP(ecModel_batch,1);
        if ~isempty(solution.x)
            load('GECKO/geckomat/parameters.mat')
            %Get U_max
            gRate = solution.x(obj);
            %Introduce curation of problematic Kcats
            ecModel_batch = curateKcatValues(ecModel_batch);
            %Fit sigma parameter and adjust protein pool upper bound
            ecModel_batch = rescaleBioMassComposition(ecModel_batch,phenotype{1},biomass_phen);
            cd ../kcat_sensitivity_analysis
            OptSigma = sigmaFitter(ecModel_batch,0.46,gRate,0.5);
            ecModel_batch.ub(end) = 0.46*OptSigma*0.5;
            %Constrain acetate production
            %ecModel_batch.ub(acEx) = 0.6;
            %Get exchange fluxes from a new solution
            solution = solveLP(ecModel_batch,1);
            GUR  = solution.x(cSource);
            EtEx = solution.x(EtOH);  
            %Get the top used enzymes
            T = topUsedEnzymes(solution.x,ecModel_batch,{'glc'},ecModelName,false);
            [limitations,breakFlag] = findTopLimitations(ecModel_batch,[]);
        else
            gRate = 0;
            GUR   = 0;
            EtEx  = 0;
        end
        cd ../utilities
        f_factor = OptSigma;
        [kcat,rxnIdx,rxnName,MW] = getKcat(ecModel_batch,T.prots_glc{1});
        modelsData = [modelsData; {modelFile phenotype{1} nRxns nMets nGnes nRxns_ec ...
                      nMets_ec nEnz_ec isoEnzymes Promiscuous Complexes ...
                      RxnWithKcat gRate GUR EtEx yeastParam.GUR yeastParam.EtOH ...
                      T.prots_glc{1} T.Usages_glc(1) kcat(1) rxnName(1) MW(1) ...
                      limitations{1} limitations{4} limitations{5} limitations{6}} f_factor];
        close all
        clc
        save(['../../../../ecModels/ec_' modelFile '_GEM/ecModel_batch.mat'],'ecModel_batch')
        %check ethanol production capabilities
        temp = setParam(ecModel_batch,'obj',EtOH,1);
        temp = setParam(temp,'lb',obj,0.1*gRate);
        solution = solveLP(temp,1);
        EtEx     = solution.x(EtOH); 
    end
end
modelsData.Properties.VariableNames ={'model' 'phenotype' 'nRxns' 'nMets' 'nGenes' 'ec_nRxns' 'ec_nMets' 'ec_nEnz' ...
            'isoenzymes' 'promiscuous' 'complexes' 'RxnWithKcat' 'gRate' ...
            'GUR' 'EtEx' 'GUR_exp' 'EtOH_exp' 'topProt' 'topUsage' 'kcat' ...
            'rxnName' 'MW' 'limEnz' 'limKcat' 'gCC' 'limRxn' 'f_factor'};
%Save results table
mkdir('../../../results')
writetable(modelsData,'../../../results/ecModels_metrics.txt','Delimiter','\t','QuoteStrings',false)
        