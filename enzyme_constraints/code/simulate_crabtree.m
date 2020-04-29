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
rxns   = {'r_1631' 'r_1549' 'r_1810' 'r_2033' 'r_1634' 'r_0659' 'r_2045_REV'}; ;'r_0662'%'r_2045_REV' 'r_0659' 
bounds = [1E-5 1E-5 1E-5 0.05 0.6 0 0];%0 0
temp   = setParam(temp,'ub',rxns,bounds);
%Set NGAM
NGAM = 0.7;
temp = setParam(temp,'lb',NGAMindex,NGAM);
temp = setParam(temp,'ub',NGAMindex,NGAM);
%Introduce manually curated Kcats
cd (current)
temp = curate_ecYeast_parameters(temp);
solution = solveLP(temp);
if ~isempty(solution.f)
    maxGrowth = solution.x(bioIndex);
    disp(['The maximum growth rate for the model on glucose minimal media is: ' num2str(maxGrowth) ' [g/gDw h]'])
    %Fit sigma factor in order to reach max experimental growth
    if maxGrowth<0.4
        cd GECKO/geckomat/kcat_sensitivity_analysis
        OptSigma = 0.48;%sigmaFitter(temp,0.46,0.4,0.5);
        temp = setParam(temp,'ub',P_Index,0.46*0.5*OptSigma);
    end
end
fermentation = false;
j=0;
i=0;
SSE      = [];
results  = [];
error    = [];
respL = GAM_0;
upper = 1.1*GAM_0;
fermL = 0.75*GAM_0;
for subopt_growth = 0:(0.4/80):0.4
    cd (current)
    tempModel = temp;   
    if ~fermentation
        GAM = upper;%abs(-respL+(respL-upper)*j/57);
        lastJ = j+1;
    else
    	GAM = abs(-upper+(upper-fermL)*(j-lastJ)/(80-lastJ));
    end
    %tempModel.S(GAMmets,bioIndex) = sign(tempModel.S(GAMmets,bioIndex))*GAM;
    %Simulate!
    cd GECKO/geckomat/utilities
    solution = simulateChemostat(tempModel,subopt_growth,[cSource bioIndex],true);
    j = j+1;
    exchangeVector = solution(exchIndexes);   
    cd(current)
    if ~fermentation
        disp(['Dilution rate = ' num2str(subopt_growth) ': Respiration'])
        if exchangeVector(4)>1E-2 
            disp(['The critical dilution rate is: ' num2str(subopt_growth) ' [1/h]'])
            fermentation = true;
        end
    else
        disp(['Dilution rate = ' num2str(subopt_growth) ': Fermentation'])
    end
    exchangeVector(exchangeVector==0) = 1E-6;
    newRow  = [subopt_growth, exchangeVector'];
    results = [results; newRow];
    
    %compare with experimental data
    %Search subopt_growth in Drate exp data
    [~,dataIndex] = ismember(subopt_growth,data.D_h_1_);
    if dataIndex>0
        Drate = data.D_h_1_;
        expExchFlux = [data.qglucose(dataIndex);data.qO2(dataIndex);data.qCO2(dataIndex);data.qethanol(dataIndex)];
        expExchFlux(expExchFlux==0) = 1E-6;
        %Get median relative error for exchange fluxes prediction for each
        %Drate
        SSE = [SSE; mean(abs(expExchFlux-exchangeVector)./(expExchFlux))];
    end
    i = i+1;
end
%Display the median relative error for exchange fluxes prediction across
%dilution rates
SSE = mean(SSE);
disp(['The median error is: ' num2str(SSE)]);
error = [error;SSE];
%Plot results
names = {'Glucose' 'Oxygen' 'CO2' 'Ethanol'};
for i=1:(length(exchIndexes))
    plot(results(:,1),results(:,i+1),'LineWidth',3)
    hold on
end
%Add experimental data points
vector = [5 3 4 6];
for i=vector
    scatter(Drate,table2array(data(:,i)))
    hold on
end
legend(names)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Exchange fluxes [mmol/gDw h]','FontSize',18)
xlim([0 max(results(:,1))])
hold off
