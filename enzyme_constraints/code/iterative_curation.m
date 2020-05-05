function model = iterative_curation(model,exp_gRate,bioY_exp,overpredicted,iterations)
if nargin<7
    iterations = 10;
end
glcIndx     = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'));
gIndx       = find(strcmpi(model.rxnNames,'biomass pseudoreaction'));
Pindex      = find(contains(model.rxnNames,'prot_pool'));
sol         = solveLP(model,1);
gSim        = sol.x(gIndx);
prediction  = (sol.x(gIndx)/(sol.x(glcIndx)*0.18));
error_gRate = (gSim-exp_gRate)/exp_gRate;
disp(['Error in growth rate prediction is: ' num2str(error_gRate)])
temp = model;
cd specific_scripts
if overpredicted
	factor    = 0.1;
	direction = -1;
	comStr    = 'limKcat = limKcat(find(coefficients>0),:);';
else
	factor    = 10;
	direction = 1;
	comStr = 'limKcat = limKcat(find(coefficients<0),:);';
end
j = 1;
modifications = {};
stopCriteria = false;
while ~stopCriteria & j<iterations
    disp(['Iteration #' num2str(j)])
    limKcat = [];
    temp = setParam(temp,'lb',gIndx,0.999999*gSim);
    temp = setParam(temp,'ub',gIndx,1.000001*gSim);
    temp = setParam(temp,'obj',gIndx,1);
    [limKcat,~] = findTopLimitationsAll(temp,[],glcIndx,factor,'ascend');
    %Filter out results   
    if ~isempty(limKcat)
        coefficients = limKcat{5};
        limKcat = table(limKcat{1},limKcat{2},limKcat{3},limKcat{4},limKcat{5},limKcat{6});
        eval(comStr)
    end
    
    if ~isempty(limKcat)
        j = j+1;
        %Modify top-ranked Kcat
        k = 1;
        stopEnz = false;
        while ~stopEnz & k<=height(limKcat)
            metIndex = table2array(limKcat(k,2));
            enzyme   = table2cell(limKcat(k,1));
            indexes  = find(temp.S(metIndex,:));
            indexes  = indexes(1:(end-1));
            if ~ismember(enzyme,modifications)
                tempCell = table2cell(limKcat(k,6));
                tempM = temp;
                tempM.S(metIndex,indexes) = temp.S(metIndex,indexes)/factor;
                sol = solveLP(tempM,1);
                if ~isempty(sol.x)
                    newPrediction = (sol.x(gIndx)/(sol.x(glcIndx)*0.18));
                    if sign(newPrediction-prediction)==sign(direction)
                        disp(enzyme{1})
                        disp(tempCell{1})
                        prediction = newPrediction;
                        disp(['The curated biomass yield is: ' num2str(prediction)])
                        modifications = [modifications; enzyme];
                        stopEnz = true;
                        temp = tempM;
                    end
                end
            end
            k = k+1;
        end
    else
        j = iterations+1;
    end
    predError = (prediction-bioY_exp)/bioY_exp;
    if overpredicted
        stopCriteria = predError <= 0.1;
    else
        stopCriteria = predError >= -0.1;
    end
end
model = temp;
end
