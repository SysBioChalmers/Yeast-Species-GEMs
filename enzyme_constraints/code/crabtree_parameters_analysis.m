%Crabtree effect parameters sensitivity analysis
%
% Kcat sensitivity analysis on ethanol and biomass yields, miu max and
% critical dilution rate for a given set of selected genes.
%
% Last modified. Ivan Domenzain 2020/04/30

clear
close all
clc
current = pwd;
%Load ecYeastGEM
load(['../ecModels/ec_Saccharomyces_cerevisiae_GEM/ecModel_batch_curated.mat'])
model = ecModel_batch;
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
%set media constraints
cd geckomat/kcat_sensitivity_analysis
temp = changeMedia_batch(model,'D-glucose exchange (reversible)');
%Block selected rxns (acetaldehyde acetate pyruvate butanediol)
rxns   = {'r_1631' 'r_1549' 'r_1810' 'r_2033' 'r_1634' 'r_0659No1' 'r_2045_REV'}; %;'r_0662'%'r_2045_REV' 'r_0659'
bounds = [1E-5 1E-5 1E-5 0.05 0.62 0 0];%0 0
temp   = setParam(temp,'ub',rxns,bounds);
cd (current)
%Set NGAM
NGAM = 0.7;
temp = setParam(temp,'lb',NGAMindex,NGAM);
model = setParam(temp,'ub',NGAMindex,NGAM);
glc_MW = 0.180156;
Et_MW  = 0.04607;
%{'FBA' 'PGK' 'PYK' 'ATP3' 'ATP5'}
%factors = [0.1 0.2 0.3 0.4 0.5 0.75 1 1.25 1.5 1.75 2 5 10 100];
factors = [1 2 3 4 5 6 7 8 9 10];
parameters = {'\mu_{max}' 'D_{crit}' 'EtOH yield' 'bio yield'};
Kcat_FC = factors';

close all
%for gene = {'FBA' 'PGK' 'PYK' 'ATP1' 'ATP3' 'ATP5'}
for gene = {'PGK'}    
    miu_max  = [];
    D_crit   = [];
    Et_yield = [];
    bioYield = [];
    for Fchange = factors
        cd GECKO/geckomat/utilities
        enzIndex = find(contains(model.enzNames,gene));
        enzyme   = model.enzymes{enzIndex};
        enzIndex = find(contains(model.metNames,enzyme));
        [kcat,rxnIdx,rxnName,MW] = getKcat(model,enzyme);
        temp = model;
        %Get model with modified enzyme activity
        temp.S(enzIndex,rxnIdx) = temp.S(enzIndex,rxnIdx)/Fchange;
        %Solve
        %Fit sigma factor in order to reach max experimental growth
        OptSigma = 0.44;%sigmaFitter(temp,0.46,0.4,0.5);
        temp = setParam(temp,'ub',P_Index,0.46*0.5*OptSigma);
        solution = solveLP(temp,1);
        if ~isempty(solution.f)
            maxGrowth = solution.x(bioIndex);
            miu_max   = [miu_max;maxGrowth];
            disp(['The maximum growth rate for the model on glucose minimal media is: ' num2str(maxGrowth) ' [g/gDw h]'])
            solution = simulateChemostat(temp,maxGrowth,[cSource bioIndex],true);
            cd (current)
            gRates = 0.5*maxGrowth:(0.4*maxGrowth/30):maxGrowth;
            [~,D_Et] = crabtree_chemostats(temp,[cSource bioIndex],exchIndexes,gRates);
            D_crit   = [D_crit;D_Et];
            Et_yield = [Et_yield; (solution(ethIndex)/solution(cSource)*(Et_MW/glc_MW))];
            bioYield = [bioYield; (solution(bioIndex)/(solution(cSource)*glc_MW))];
        else
            disp(['Not feasible ' gene{1} ' ' num2str(Fchange)])
            miu_max = [miu_max;0];
            D_crit   = [D_crit;0];
            Et_yield = [Et_yield;0];
            bioYield = [bioYield;0];
        end
        clc
        cd (current)
    end
    eval([gene{1} '=table(Kcat_FC,miu_max,D_crit,Et_yield,bioYield);'])
    eval(['geneResults = table2array(' gene{1} ');'])
    for i=1:5
        geneResults(:,i) = (geneResults(:,i))/geneResults(1,i);
    end
    figure
    plot(geneResults(:,1),geneResults(:,2),geneResults(:,1),geneResults(:,3),...
         geneResults(:,1),geneResults(:,4),geneResults(:,1),geneResults(:,5),...
         'LineWidth',4)
    xlim([1 10])
    ylim([0.9 1.1])
    %set(gca,'XScale', 'log','FontSize',22)
    %legend(parameters,'FontSize',18)
    title(gene{1},'FontSize',18)
    saveas(gcf,['../results/Crabtree_Kcat_sensitivity/crabtree_sensitivity_' gene{1} '.jpg'])
    savefig(['../results/Crabtree_Kcat_sensitivity/crabtree_sensitivity_' gene{1} '.fig'])
end    