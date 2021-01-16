%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OptSigma = sigmaFitter(model_batch,Ptot,expVal,f)
% 
% Function that fits the average enzyme saturation factor in an ecModel
% according to a provided experimentally measured value for the objective
% function (i.e. growth rate at specified conditions)
%
% INPUTS:
%       model_batch     An EC batch model with an initial sigma factor
%                       assigned
%       Ptot            Total protein amount in the model (Experimental)
%                       [g/gDw]
%       expVal          Experimentally measured value for the objective function
%       f               Estimated mass fraction of enzymes in model [g/g]
%
% OUTPUTS:
%       optSigma    The optimal sigma value obtained
%
% Benjamin Sanchez      2018-08-10
% Ivan Domenzain        2020-04-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptSigma,error] = sigmaFitter(model_batch,Ptot,expVal,f,objPos)
disp('Fitting protein pool parameters')
objValues   = [];
errors    = [];
sigParam  = [];
poolIndex = find(strcmpi(model_batch.rxnNames,'prot_pool_exchange'));
if nargin<5
    objPos = find(model_batch.c);
end
%Relax bounds for the objective function
model_batch.lb(objPos) = 0;
model_batch.ub(objPos) = 1000;

for i=1:100
    %Constrains the ecModel with the i-th sigma factor
    sigma = i/100;
    model_batch.ub(poolIndex) = sigma*Ptot*f; 
    solution  = solveLP(model_batch,1);
    if isempty(solution.x)
        solution.x=zeros(length(model_batch.rxns),1);
    end
    objValues = [objValues; solution.x(objPos)];
    error     = abs(((expVal-solution.x(objPos))/expVal)*100);
    errors    = [errors; error];
    error     = num2str(((expVal-solution.x(objPos))/expVal)*100);
    sigParam  = [sigParam; sigma];
    %disp(['Fitting sigma factor: ' num2str(sigma) '   error: ' error '%'])
end
if min(errors)==100
    indexes = find(errors~=100);
    minIndx = min(indexes);
else
    [~, minIndx] = min(errors);
end
OptSigma  = sigParam(minIndx);
error     = errors(minIndx);
disp(['Fitted sigma factor: ' num2str(OptSigma) '   error: ' num2str(errors(minIndx)) '%'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
