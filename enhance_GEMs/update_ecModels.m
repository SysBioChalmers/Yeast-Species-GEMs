%update_ecModels
%
% Ivan Domenzain. 2020-03-04
%

%Clone GECKO repository
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('pull')
%Locate the correct branch
git('checkout fix/updateDatabases') 
clc
fileNames = dir('../models');
%Upload KEGG codes info
fID        = fopen('../../ComplementaryData/yeasts_keggCodes.txt');
keggCodes  = textscan(fID,'%s %s %s','Delimiter','\t','HeaderLines',1);
variables  = {'strain','modelName','KEGG'};
keggCodes = table(keggCodes{1},keggCodes{2},keggCodes{3},'VariableNames', variables);
%Substitute some scripts in GECKO
copyfile('../updateDatabases.m','geckomat/get_enzyme_data/updateDatabases.m')
current = pwd;
for i=1:length(fileNames)
    cd(current)
    file = fileNames(i).name;
    if contains(file,'.mat')
        disp(file)
        modelName = file(1:(end-4));
        disp(modelName)
        DBname    = ['../../ComplementaryData/databases/uniprot/' modelName '.tab'];
        %If a uniprot file is available for this organism in the repository
        %then paste it in GECKO
        if isfile(DBname)
            copyfile(DBname,'databases/uniprot.tab')
            index = find(strcmpi(keggCodes.modelName,modelName),1);
            keggCode = keggCodes.KEGG(index);
            disp(keggCode)
            %Generate protDatabase (uniprot and KEGG for the desired
            %organism)
            cd geckomat/get_enzyme_data
            updateDatabases('');
            copyfile('../../Databases/protDatabase.mat',['../../Databases/' modelName '_protDatabase.mat'])
            load(['../../../models/' file])
            cd ..
            ecModelName = ['ec_' modelName '_GEM'];
            [ecModel,ecModel_batch] = enhanceGEM(reducedModel,'COBRA','ecRhtoGEM','v.1.0');
        end
    end
end