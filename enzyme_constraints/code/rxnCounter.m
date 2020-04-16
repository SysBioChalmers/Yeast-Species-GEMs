%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [isoEnzymes,Promiscuous,Complexes,RxnWithKcat] =  rxnCounter(model,originalModel)
%
% Function that receives an EC model and counts the number of reactions 
% with isoenzymes, enzyme complexes, Promiscuous enzymes and rxns with 
% Kcat values.
%
%Ivan Domenzain.       Last edited: 2018-08-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isoEnzymes,Promiscuous,Complexes,RxnWithKcat] =  rxnCounter(model,originalModel)
    % Get enzymes indxs
    origRxns    = originalModel.rxns;
    enzMetIndxs = find(contains(model.metNames,'prot_'));
    enzNames    = model.metNames(find(contains(model.metNames,'prot_')));
    enzRxnIndxs = find(contains(model.rxnNames,'prot_'));
    enzRxnIndxs = enzRxnIndxs(3:end);
    isoEnzymes  = find(contains(model.rxnNames,'(arm)'));
    isoEnzymes  = length(isoEnzymes);
    S_enzymes   = model.S(enzMetIndxs,1:enzRxnIndxs(1)-1);
    [m, n]      = size(S_enzymes);
    promDist    = [];
    promNames   = [];
    compDist    = [];
    isoEDist    = [];
    ylab        = 'Frequency';
    
    % Isoenzymes 
    for i=1:length(origRxns)
        ID = origRxns{i};
        ID = [ID 'No'];
        nIso = sum(contains(model.rxns,ID));
        if nIso>1
            isoEDist = [isoEDist; nIso];
        end  
    end
    nBins = 20;
    titlestr = [num2str(isoEnzymes) ' Reactions with isoenzymes'];
    plotHist(isoEDist,'# of isoenzymes per reaction',ylab,titlestr,nBins);
    % Promiscuous 
    Promiscuous = 0;
    for i=1:m
        nProm = length(find(S_enzymes(i,:)~=0));
        %Indx = enzMetIndxs(i);
        if nProm>1
            Promiscuous = Promiscuous+1;
            promNames   = [promNames;enzNames(i)];
            promDist    = [promDist; nProm];
        end       
    end
    nBins = 50;
    titlestr = [num2str(Promiscuous) ' Promiscuous enzymes'];
    [~,I]    = max(promDist);
    promDist(I) = [];
    [~,I]    = max(promDist);
    promDist(I) = [];
    disp(['The most promiscuous enzyme is: ' promNames{I}])
    %plotHist(promDist,'# of reactions per enzyme',ylab,titlestr,nBins);
    
    % Complexes 
    Complexes = 0; RxnWithKcat = 0;
    for i=1:n
        nComp = length(find(S_enzymes(:,i)~=0));
        if nComp>0
            RxnWithKcat = RxnWithKcat+1;
            if nComp>1
                Complexes = Complexes+1;
                compDist  = [compDist; nComp]; 
            end
        end 
    end
    nBins = 20;
    titlestr = [num2str(Complexes) ' Enzyme complexes'];
    %plotHist(compDist,'# of subunits per complex',ylab,titlestr,nBins);    
end

function plotHist(var,xlab,ylab,titlestr,nBins)
    figure
    hist(var,nBins)
    title(titlestr,'FontSize',18)
    xlabel(xlab,'FontSize',18)
    ylabel(ylab,'FontSize',18)
end

