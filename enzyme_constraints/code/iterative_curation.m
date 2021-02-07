function [model,prediction] = iterative_curation(model,expVal,underpredicted,curIndx,fixIndx,objIndx,factor,iterations)
if nargin<7
	iterations = 10;
end
sol         = solveLP(model,1);
fixVal      = sol.x(fixIndx);
prediction  = sol.x(curIndx);
%error_gRate = (simVal-expVal)/expVal;
%disp(['Error in growth rate prediction is: ' num2str(error_gRate)])
temp = model;
%Fix growth rate value and objective coefficient
temp = setParam(temp,'lb',fixIndx,0.99*fixVal);
temp = setParam(temp,'ub',fixIndx,1.01*fixVal);
temp = setParam(temp,'obj',objIndx,1);
cd specific_scripts
if underpredicted
	%factor    = 0.1;
	direction = 1;
	comStr    = 'limKcat = limKcat(find(coefficients>0),:);';
    sorting   = 'descend';
else
	%factor    = 10;
	direction = -1;
	comStr = 'limKcat = limKcat(find(coefficients<0),:);';
    sorting   = 'ascend';
end
j = 1;
modifications = {};
stopCriteria  = false;
noLims        = false;
while ~stopCriteria & j<iterations & ~noLims
    disp(['Iteration #' num2str(j)])
    limKcat = [];
    %Calculate ECC
    [limKcat,noLims] = findTopLimitationsAll(temp,[],curIndx,factor,sorting);
    %Filter out results (negative or positive)
    if ~isempty(limKcat)
        coefficients = limKcat.ECC;
        %limKcat = table(limKcat{1},limKcat{2},limKcat{3},limKcat{4},limKcat{5},limKcat{6});
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
                    newPrediction = (sol.x(curIndx));
                    if sign(newPrediction-prediction)==sign(direction)
                        disp(enzyme{1})
                        disp(tempCell{1})
                        prediction = newPrediction;
                        disp(['The curated value is: ' num2str(prediction)])
                        modifications = [modifications; enzyme];
                        stopEnz = true;
                        temp = tempM;
                    end
                end
            end
            if k==height(limKcat)
                noLims = true;
            end
            k = k+1;
        end
    else
        j = iterations+1;
    end
    predError = (prediction-expVal)/expVal;
    if underpredicted
        stopCriteria = predError >= -0.1;     
    else
        stopCriteria = predError <= 0.1;
    end
end
model = temp;
end
