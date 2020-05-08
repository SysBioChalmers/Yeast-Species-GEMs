function model = iterative_curation(model,expVal,expYield,overpredicted,in_idx,out_idx,iterations,MW)
if nargin<8
    MW = 1;
    if nargin<7
        iterations = 10;
    end
end
sol         = solveLP(model,1);
simVal      = sol.x(out_idx);
prediction  = (sol.x(out_idx)/(sol.x(in_idx)*MW));
error_gRate = (simVal-expVal)/expVal;
disp(['Error in growth rate prediction is: ' num2str(error_gRate)])
temp = model;
%Fix growth rate value and objective coefficient
temp = setParam(temp,'lb',out_idx,0.999999*simVal);
temp = setParam(temp,'ub',out_idx,1.000001*simVal);
temp = setParam(temp,'obj',out_idx,1);
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
    %Calculate ECC
    [limKcat,~] = findTopLimitationsAll(temp,[],in_idx,factor,'ascend');
    %Filter out results (negative or positive)
    if ~isempty(limKcat)
        coefficients = limKcat{5};
        limKcat = table(limKcat{1},limKcat{2},limKcat{3},limKcat{4},limKcat{5},limKcat{6});
        eval(comStr)
    end
    %Modify top-ranked Kcat
    if ~isempty(limKcat)
        j = j+1;
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
                    newPrediction = (sol.x(out_idx)/(sol.x(in_idx)*MW));
                    if sign(newPrediction-prediction)==sign(direction)
                        disp(enzyme{1})
                        disp(tempCell{1})
                        prediction = newPrediction;
                        disp(['The curated yield is: ' num2str(prediction)])
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
    predError = (prediction-expYield)/expYield;
    if overpredicted
        stopCriteria = predError <= 0.1;
    else
        stopCriteria = predError >= -0.1;
    end
end
model = temp;
end
