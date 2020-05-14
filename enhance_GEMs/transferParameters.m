function transferParameters(parameters,model,modelName)
orgNames   = parameters{2};
modelNames = parameters{3};
keggCodes  = parameters{4};
gRates     = parameters{5};
cd GECKO/geckomat/limit_proteins
Ptot     = sumProtein(model);
index    = find(strcmpi(modelNames,modelName),1);
keggCode = keggCodes{index};
org_name = lower(orgNames{index});
gRate    = gRates(index);
yeastParam = struct('Ptot',Ptot,'gR_exp',gRate,'org_name',org_name,'keggID',keggCode);
save('../parameters.mat','yeastParam')
cd ../..
end