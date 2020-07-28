%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveSSModel(model,upDATE)
% Saves model as a .xml, .txt and .yml file. Also updates complementary
% files (boundaryMets.txt, README.md and dependencies.txt).
%
% model     model structure to save (note: must be in COBRA format)
% upDATE    logical =true if updating the date in the README file is needed
%           (opt, default true)
%
% Benjamin J. Sanchez
% Feiran Li -2019.11.12 - Modify to fit Specific models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveSSModel(model,upDATE)

if nargin < 2
    upDATE = true;
end

%Get and change to the script folder, as all folders are relative to this
%folder
scriptFolder = fileparts(which(mfilename));
currentDir = cd(scriptFolder);

%Set minimal media
cd Reconstruction/otherchanges/
model = minimal_Y6(model);


%Delete model.grRules (redundant and possibly conflicting with model.rules):
 if isfield(model,'grRules')
     model = rmfield(model,'grRules');
 end

%Update SBO terms in model:
% model = addSBOterms(model);
cd ../../
%Check if model is a valid SBML structure:
writeCbModel(model,'sbml','tempModel.xml');
 [~,errors] = TranslateSBML('tempModel.xml');
if ~isempty(errors)
    delete('tempModel.xml');
    error('Model should be a valid SBML structure. Please fix all errors before saving.')
end

%Update .xml, .txt and .yml models:
modelName = split(model.id,' specific');
modelName = modelName{1};
copyfile('tempModel.xml',['ModelFiles/xml/',modelName,'.xml'])
delete('tempModel.xml');
writeCbModel(model,'text',['ModelFiles/txt/',modelName,'.txt']);
%exportForGit(model,modelName,'..',{'yml'});

%Update README file: date + size of model
copyfile('../README.md','backup.md')
fin  = fopen('backup.md','r');
fout = fopen('../README.md','w');
still_reading = true;
while still_reading
    inline = fgets(fin);
    if ~ischar(inline)
        still_reading = false;
    else
        if startsWith(inline,'* Last update: ') && upDATE
            inline = ['* Last update: ' datestr(datetime,'yyyy-mm-dd') newline];
        elseif startsWith(inline,['|_',modelName,'_|'])
            inline = [['|_',modelName,'_|Yeast8.3|']...
                num2str(length(model.rxns)) '|' ...
                num2str(length(model.mets)) '|' ...
                num2str(length(model.genes)) '|'...
                datestr(datetime,'yy-mm-dd') '|'...
                newline];
        end
        fwrite(fout,inline);
    end
end
fclose('all');
delete('backup.md');

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
copyfile(['ModelFiles/xml/',modelName,'.xml'],'backup.xml')
fin  = fopen('backup.xml','r');
fout = fopen(['ModelFiles/xml/',modelName,'.xml'],'w');
still_reading = true;
while still_reading
    inline = fgets(fin);
    if ~ischar(inline)
        still_reading = false;
    else
        if ~isempty(regexp(inline,'[0-9]e-?00[0-9]','once'))
            inline = regexprep(inline,'(?<=[0-9]e-?)00(?=[0-9])','0');
        end
        fwrite(fout,inline);
    end
end
fclose('all');
delete('backup.xml');

%Switch back to original folder
cd(currentDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
