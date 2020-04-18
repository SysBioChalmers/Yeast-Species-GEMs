%FSEOF on ethanol production at high growth rates
current = pwd;
%Clone GECKO repository
%git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
%git('pull')
%Locate the correct branch
git('stash')
git('checkout feat/add_FSEOF_utilities') 
clc
cd ..
%Get relevant rxns ids and indexes
CS_MW    = 0.18015;
c_source = 'D-glucose exchange (reversible)';

%Retrieve ecModel names
fileNames = dir('../ecModels');
mkdir('../results/EtOH_production')
for i=1:length(fileNames)
    cd(current)
    file = fileNames(i).name;
    if contains(file,'ec_')
        name = strrep(file,'ec_','');
        name = strrep(name,'_GEM','');
        disp(file)
        load(['../ecModels/' file '/ecModel_batch_curated.mat'])
        model  = ecModel_batch;
        target = model.rxns(find(strcmpi(model.rxnNames,'ethanol exchange')));
        growth_pos = find(model.c);
        CS_index   = find(strcmpi(model.rxnNames,c_source));
        WT_solution  = solveLP(model,1);
        WT_MAX_yield = WT_solution.x(growth_pos)/(WT_solution.x(CS_index)*CS_MW);
        disp(['The biomass yield is ' num2str(WT_MAX_yield) '[g biomass/g carbon source]'])
        cd GECKO/geckomat/utilities/ecFSEOF
        output  = ['../../../../../results/EtOH_production/' name '_EtOH_genes.txt'];
        results = run_ecFSEOF(model,target,c_source,[0.99*WT_MAX_yield WT_MAX_yield],16,output);
    end
end