%update_ecModels
%
% Ivan Domenzain. 2020-03-11
%
current = pwd;
%Clone GECKO repository
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git('pull')
%Locate the correct branch
git('checkout fix/updateDatabases') 
clc
fileNames = dir('../../models');
%Upload organism and model specific parameters
fID         = fopen('../../data/yeasts_parameters.txt');
yeastsParam = textscan(fID,'%s %s %s %s %f %f %f %s','Delimiter','\t','HeaderLines',1);
%Replace scripts in GECKO:
scripts = dir('../specific_scripts');
for i = 1:length(scripts)
    script = scripts(i).name;
    if contains(script,'.m')
        fullName   = ['../specific_scripts/' script];
        %Retrieve script path within GECKO
        GECKO_path = dir(['**/' script]);
        if ~isempty(GECKO_path)
            GECKO_path = GECKO_path.folder;
            %Replace script in GECKO in its container subfolder
            copyfile(fullName,GECKO_path)
        end
    end
end
delete databases/chemostatData.tsv
delete databases/prot_abundance.txt
mkdir('../../ecModels')
for i=1:length(fileNames)
    cd(current)
    file = fileNames(i).name;
    if contains(file,'.mat')
        modelName = file(1:(end-4));
        disp(modelName)
        %Load model
        load(['../models/' file])
        %Convert to RAVEN format
        if isfield(reducedModel,'rules')
            model = ravenCobraWrapper(reducedModel);
        else
            model = reducedModel;
        end
        %Get oxPhos GPRs from yeast GEM
        cd specific_scripts
        model = getOxPhosGPRs(model);
        cd ..
        %Transfer model parameters to GECKO
        transferParameters(yeastsParam,model,modelName)       
        %If a uniprot file is available for this organism in the repository
        %then paste it in GECKO
        DBname = ['../../../ComplementaryData/databases/uniprot/' modelName '.tab'];
        if isfile(DBname)
            copyfile(DBname,'databases/uniprot.tab')
            %Generate protDatabase (uniprot and KEGG for the desired
            %organism)
            cd geckomat/get_enzyme_data
            updateDatabases('');
            copyfile('../../Databases/protDatabase.mat',['../../Databases/' modelName '_protDatabase.mat'])
            ecModelName = ['ec_' modelName '_GEM'];        
            mkdir(['../../models/' ecModelName])
            cd ..
            [ecModel,ecModel_batch] = enhanceGEM(model,'COBRA',ecModelName,'v.1.0');
            save(['../models/' ecModelName '/ecModel.mat'],'ecModel')
            save(['../models/' ecModelName '/ecModel_batch.mat'],'ecModel_batch')
            newDir = ['../../../ecModels/' ecModelName];
            mkdir(newDir)
            movefile(['../models/' ecModelName '/*'],newDir)
        end
    end
end