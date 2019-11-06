findCorebiomass()
for i = 1:length(StrianData.genes)
    if find(contains(model.genes,StrianData.genes{i}))
    model_temp  = removeGenes(model,StrianData.genes{i},true,true,true);
    
    bio_rxn = {'r_4048';'r_4049';'r_4050';'r_4598';'r_4599'};% all biomass pseudoreactions except protein, we will manually add amino acid production, since the precursor in protein_pseudoreaction are aa_chargerd tRNAs
    [~,bio_rxn_index] = ismember(bio_rxn,model_temp.rxns);
    mets = [];
    for j = 1:length(bio_rxn_index)
        mets_temp = find(model_temp.S(:,bio_rxn_index(j))< 0 );
        mets = [mets;mets_temp];
    end
    
    results = canProduce(model_temp,model_temp.mets(mets));
    if all(results)
        %model = model_temp;
    else
        for j = 1:length(results)
            if ~logical(results(j))
                j
            reason = [reason;StrianData.genes{i},mets_temp.mets(mets(j))];
            end
        end
    end
    end
end


    
