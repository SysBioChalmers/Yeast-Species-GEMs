function model = addrxnBack(model,model_original,rxnsList,GPR)
if isfield(model,'rxnMiriams')
    model = ravenCobraWrapper(model);
end
if nargin < 4
    GPR = cell(length(rxnsList),1);
end
for j = 1:length(rxnsList)
    [~,ID] = ismember(rxnsList(j),model_original.rxns);
    [~,tmp] = ismember(rxnsList(j),model.rxns);
    if isempty(tmp)
        warning(['model alrady have this rxn:',rxnsList(j)])
    end
    rxn.rev =  model_original.lb(ID);
    if rxn.rev == -1000
        rxn.rev = 1;
    else
        rxn.rev = 0;
    end
    mets = find(model_original.S(:,ID));
    matrix.coef = model_original.S(mets,ID);
    matrix.metID = model_original.mets(mets);
    matrix.metNames = model_original.metNames(mets);
    matrix.metKEGGID = model_original.metKEGGID(mets);
    matrix.metChEBIID = model_original.metChEBIID(mets);    
    matrix.metFormulas = model_original.metFormulas(mets);
    matrix.metCharges = model_original.metCharges(mets);
    matrix.metFormulas = model_original.metFormulas(mets);
    matrix.metCharges = model_original.metCharges(mets);
    
    %mapping mets to model.metnames, get s_ index for new mets
    [~,metindex] = ismember(matrix.metNames,model.metNames);
    if any(metindex == 0)
        for k = 1:length(metindex)
            if metindex(k) == 0
                model = addMetabolite(model,matrix.metID{k}, ...
                    'metName',matrix.metNames{k});
            end
        end
    end
    
for i = 1:length(matrix.metNames)
    [~,metID] = ismember(matrix.metNames(i),model.metNames);
    if metID ~= 0
        model.metFormulas{metID} = matrix.metFormulas{i};
        model.metCharges(metID)  =  matrix.metCharges(i);
        model.metKEGGID{metID}   = matrix.metKEGGID{i};
        model.metChEBIID{metID}  = matrix.metChEBIID{i};
        model.metNotes{metID}    = 'NOTES: added after Gapfilling; ';
    end
end
    %add reaction back
    rxnID   = model_original.rxns{ID};
    [model,rxnIndex] = addReaction(model,...
        rxnID,...
        'reactionName', rxnID,...
        'metaboliteList',matrix.metID,...
        'stoichCoeffList',matrix.coef,...
        'reversible',rxn.rev,...
        'geneRule',GPR{j},...
        'checkDuplicate',1);
    if isempty(rxnIndex)
        [~,rxnIndex] = ismember(rxnID,model.rxns);
    end
    % Add rxn annotation:
    model.rxnNames{rxnIndex}      = model_original.rxnNames{ID};
    model.rxnECNumbers(rxnIndex)  = model_original.rxnECNumbers(ID);
    model.rxnKEGGID(rxnIndex)     = model_original.rxnKEGGID(ID);
    model.rxnMetaNetXID(rxnIndex) = model_original.rxnMetaNetXID(ID);
    model.rxnConfidenceScores(rxnIndex) = 0;
    model.rxnNotes{rxnIndex} = 'NOTES: added after Gapfilling; ';
    % Display the 
    disp(['rxn: ',rxnID, ' has been added back to the model '])
    printRxnFormula(model,'rxnAbbrList',rxnID,'metNameFlag',true);
end

% add gene standard name for new genes
fid = fopen('../../data/databases/SGDgeneNames.tsv');
yeast_gene_annotation = textscan(fid,'%s %s','Delimiter','\t','HeaderLines',1);
fclose(fid);

geneIndex = zeros(1,1);
for i = 1: length(model.genes)
    geneIndex = strcmp(yeast_gene_annotation{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(yeast_gene_annotation{2}{geneIndex})
        model.geneNames{i} = yeast_gene_annotation{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end
model.proteins(length(model.genes)) = {''};
end