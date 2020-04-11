function transferParameters(parameters,model,modelName)
orgNames   = parameters{2};
modelNames = parameters{3};
keggCodes  = parameters{4};
gRates     = parameters{5};
GURs       = parameters{6};
EtOH_vExs  = parameters{7};
cd GECKO/geckomat/limit_proteins
Ptot     = sumProtein(model);
indexes  = find(strcmpi(modelNames,modelName),1);
keggCode = keggCodes{indexes(1)};
org_name = lower(orgNames{indexes(1)});
gRate    = mean(gRates(indexes));
GUR      = mean(GURs(indexes));
EtOH_vEx = mean(EtOH_vExs(indexes));
yeastParam = struct('Ptot',Ptot,'gR_exp',gRate,'org_name',org_name,...
                    'keggID',keggCode,'GUR',GUR,'EtOH',EtOH_vEx);
save('../parameters.mat','yeastParam')
cd ../..
end