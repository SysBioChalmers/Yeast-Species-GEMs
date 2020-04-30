%Simulate crabtree effect
clear
close all
clc
current = pwd;
%Load ecYeastGEM
load(['../ecModels/ec_Saccharomyces_cerevisiae_GEM/ecModel_batch_curated.mat'])
model = ecModel_batch;
%Load chemostat data
%Load experimental data
fileName = ('../data/chemostatData.txt');
data = readtable(fileName);

%Clone GECKO and go to utilities branch
%git('clone https://github.com/SysBioChalmers/GECKO')
cd GECKO
%git('stash')
%git('checkout feat-add_utilities')
%Get relevant indexes
bioIndex    = find(strcmpi(model.rxnNames,'biomass pseudoreaction'));
cSource     = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'));
oxyIndex    = find(strcmpi(model.rxnNames,'oxygen exchange (reversible)'));
CO2Index    = find(strcmpi(model.rxnNames,'carbon dioxide exchange'));
ethIndex    = find(strcmpi(model.rxnNames,'ethanol exchange'));
NGAMindex   = find(contains(model.rxnNames,'maintenance'));
P_Index     = find(strcmpi(model.rxnNames,'prot_pool_exchange'));
exchIndexes = [cSource;oxyIndex;CO2Index;ethIndex];
%Show GAM
%Get GAM related indexes
GAMcompounds = {'ATP' 'phosphate' 'ADP' 'H+' 'H2O'};
bioMets = find(model.S(:,bioIndex));
[~,iB]  = ismember(GAMcompounds,model.metNames(bioMets));
GAMmets = bioMets(iB);
GAM_0 = abs(model.S(GAMmets(1),bioIndex));
%set media constraints
cd geckomat/kcat_sensitivity_analysis
temp = changeMedia_batch(model,'D-glucose exchange (reversible)');
%Block selected rxns (acetaldehyde acetate pyruvate butanediol)
rxns   = {'r_1631' 'r_1549' 'r_1810' 'r_2033' 'r_1634' 'r_0659No1' 'r_2045_REV'}; %;'r_0662'%'r_2045_REV' 'r_0659'
bounds = [1E-5 1E-5 1E-5 0.05 0.62 0 0];%0 0
temp   = setParam(temp,'ub',rxns,bounds);
%Set NGAM
NGAM = 0.7;
temp = setParam(temp,'lb',NGAMindex,NGAM);
model = setParam(temp,'ub',NGAMindex,NGAM);
%Introduce manually curated Kcats
%temp = curate_ecYeast_parameters(temp);

genes = {'FBA' 'PGK' 'PYK' 'ATP3' 'ATP5'};
for gene = genes
    for factor = [0.1 0.5 1 2 5]
        cd (current)
        cd GECKO/geckomat/utilities
        enzIndex = find(contains(model.enzNames,gene));
        enzyme   = model.enzymes{enzIndex};
        enzIndex = find(contains(model.metNames,enzyme));
        [kcat,rxnIdx,rxnName,MW] = getKcat(model,enzyme);
        temp = model;
        temp.S(enzIndex,rxnIdx) = temp.S(enzIndex,rxnIdx)/factor;
        %Solve
        solution = solveLP(temp);
        if ~isempty(solution.f)
            maxGrowth = solution.x(bioIndex);
            disp(['The maximum growth rate for the model on glucose minimal media is: ' num2str(maxGrowth) ' [g/gDw h]'])
            %Fit sigma factor in order to reach max experimental growth
            if maxGrowth<0.4
                cd ../kcat_sensitivity_analysis
                OptSigma = 0.44;%sigmaFitter(temp,0.46,0.4,0.5);
                temp = setParam(temp,'ub',P_Index,0.46*0.5*OptSigma);
                cd ../utilities
            end
        else
            disp(['Not feasible ' gene{1} ' ' num2str(factor)])
        end
        gRates = 0:(.36/40):0.36;
        cd(current)
        results = crabtree_chemostats(model,[cSource bioIndex],exchIndexes,gRates);
        %Plot results
        figure
        names = {'Glucose' 'Oxygen' 'CO2' 'Ethanol'};
        for i=1:(length(exchIndexes))
            plot(results(:,1),results(:,i+1),'LineWidth',3)
            hold on
        end
        if factor>1E-6
            str = num2str(factor);
        else
            str = '0';
        end
        xlim([0 0.35])
        ylim([0 25])
        set(gca,'FontSize',22)
        title([gene{1} ': ' str '*Kcat'],'FontSize',18)
        %xlabel('Dilution rate [1/h]','FontSize',18)
        %ylabel('Exchange fluxes [mmol/gDw h]','FontSize',18)
        %legend(names,'FontSize',18)
        hold off
        cd (current)
        saveas(gcf,['../results/Crabtree_Kcat_sensitivity/crabtree_' gene{1} '_' str '.jpg'])
    end
end

