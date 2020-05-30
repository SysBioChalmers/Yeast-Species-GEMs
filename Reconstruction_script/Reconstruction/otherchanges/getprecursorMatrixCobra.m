function [proMarix,rxnMatrix,mets_test,cannotProduce,solresult] = getprecursorMatrixCobra(model_original,strains,filefolder,met,checkmode)
% This function is to get the rxnMatrix for rxn existence in the panmodel;
% and also generate the matrix for biomass precursors production.
% panmodel is raven format
% met is optional input if no met then use all mets_test
% Feiran Li 2019.12.11

% Find all components of biomass, one thing to be noted here is that
% protein is made of aa-tRNA so we need to use a seperate mets for them.


if ~exist('checkmode','var') || isempty(checkmode)
    checkmode = 1; % which will test whether the model can produce met input; checkmode == 0  means no precursor will be test, will only generate the rxnMatrix
end

if ~exist('met','var') || isempty(met)
    bio_rxn = {'r_4048';'r_4049';'r_4050';'r_4598';'r_4599';'r_4063';'r_4065'};% all biomass pseudoreactions except protein, we will manually add amino acid production, since the precursor in protein_pseudoreaction are aa_chargerd tRNAs
    [~,bio_rxn_index] = ismember(bio_rxn,model_original.rxns);
    mets_test = [];
    for j = 1:length(bio_rxn_index)
        mets_temp = find(model_original.S(:,bio_rxn_index(j))< 0 );
        mets_test = [mets_test;model_original.mets(mets_temp)];
    end
    %aa = {'s_1267','s_0956','s_0966','s_0970','s_0974','s_0982','s_1000','s_0992','s_1004','s_1007','s_1017','s_1022','s_1026','s_1030','s_1033','s_1036','s_1041','s_1046','s_1049','s_1052','s_1057'};
    aa = {'s_1267[e]','s_0956[e]','s_0966[e]','s_0970[e]','s_0974[e]','s_0982[e]','s_1000[e]','s_0992[e]','s_1004[e]','s_1007[e]','s_1017[e]','s_1022[e]','s_1026[e]','s_1030[e]','s_1033[e]','s_1036[e]','s_1041[e]','s_1046[e]','s_1049[e]','s_1052[e]','s_1057[e]'};
    mets_test = [mets_test;transpose(aa)];
else
    mets_test = met;
end

proMarix = zeros(length(strains),length(mets_test));
solresult = zeros(length(strains),length(mets_test));
rxnMatrix = zeros(length(strains),length(model_original.rxns));
path = pwd;

for i = 1:length(strains)
    fprintf([strains{i},' : No.',num2str(i),'\n']);
    cd(filefolder)
    load([strains{i},'.mat'])
    model = reducedModel;
    if isfield(model,'csense')
        model = rmfield(model,'csense');
    end

    sol = solveLP(model);
    cd(path)
    if ~isempty(sol.x) && checkmode
        [~, presentMets,solTotal] = PrecursorCheck(model,false,false,[],mets_test);
        [results,~] = ismember(mets_test,presentMets);
        proMarix(i,:) = transpose(results);% strains x mets
        solresult(i,:) = transpose(solTotal);% strains x mets
    end
    [index,~] = ismember(model_original.rxns,model.rxns);
    rxnMatrix(i,:) = transpose(index);% strains x rxn
end

cannotProduce = [];
if checkmode
    for j = 1:length(mets_test)
        strianlist = strains(proMarix(:,j) == 0);
        rep = repmat(mets_test(j),length(strianlist),1);
        cannotProduce = [cannotProduce;[strianlist,rep]];
    end
end

end
