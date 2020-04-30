function [results,D_crit] = crabtree_chemostats(model,objIndexes,exchIndexes,gRates)
cd GECKO/geckomat/utilities
fermentation = false;
results      = [];
D_crit       = [];
for subopt_growth = gRates
    tempModel = model;
    %Simulate!
    solution = simulateChemostat(tempModel,subopt_growth,objIndexes,true);
    exchangeVector = solution(exchIndexes);
    if ~fermentation
        disp(['Dilution rate = ' num2str(subopt_growth) ': Respiration'])
        if exchangeVector(4)>1E-2
            disp(['The critical dilution rate is: ' num2str(subopt_growth) ' [1/h]'])
            fermentation = true;
        end
    else
        disp(['Dilution rate = ' num2str(subopt_growth) ': Fermentation'])
        if isempty(D_crit)
            D_crit = subopt_growth;
        end
    end
    exchangeVector(exchangeVector==0) = 1E-6;
    newRow  = [subopt_growth, exchangeVector'];
    results = [results; newRow];
end
if isempty(D_crit)
    D_crit = 0;
end
end