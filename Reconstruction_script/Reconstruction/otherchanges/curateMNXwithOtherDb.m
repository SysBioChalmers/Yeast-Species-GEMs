function [MNXrxn,MNXmet] = curateMNXwithOtherDb
% This one generate the MNXrxn database with rev from seed and metacyc get
% rxnNames from kegg and get chebi IDs KEGG ids for MNXmet
% reactions.tsv from modelSEED database is downloaded directly from GitHub repository ModelSEEDDatabase
% link: github.com/ModelSEED/ModelSEEDDatabase
% reac_prop.tsv is downloaded from METAnetX database
% MNXref.mat is loaded from RAVEN toolobx feature branch 'add_MetaNetX'
% 'metaCycRxns.mat' and 'keggRxns.mat' are loaded from RAVEN/external

%% rxn rev part
current_path = pwd;
downloadMNXdb('reac_prop',current_path) %MNX reactions databse with equation
format = repmat('%s ',1,6);
format = strtrim(format);
rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);

for i = 1:length(rxn_temp)
    MNXrxn(:,i) = rxn_temp{i};
end
commentLines = startsWith(MNXrxn(:,1),'#');
MNXrxn(commentLines,:) = []; % MNX_ID	Equation	Description	Balance	EC	Source
MNXrxn(:,2) = strrep(MNXrxn(:,2),'=','<=>');
% load MNX databse ref ID mapping
% find RAVEN toolbox
IDfile = 'ravenCobraWrapper.m';
%Try to find root of toolbox:
try
    toolboxPath = which(IDfile);                %full file path
    slashPos    = getSlashPos(toolboxPath);
    toolboxPath = toolboxPath(1:slashPos(end)); %folder path
    %Go up until the root is found:
    D = dir(toolboxPath);
    while ~ismember({'.git'},{D.name})
        slashPos    = getSlashPos(toolboxPath);
        toolboxPath = toolboxPath(1:slashPos(end-1));
        D = dir(toolboxPath);
    end
    cd(toolboxPath);
catch
    disp([toolbox ' toolbox cannot be found'])
end
load('external/MetaNetX/MNXref.mat'); % requires RAVEN to checkout at feature branch of 'add_Metanetx'
% load Metacyc rxns for rev
metacycRxn = load('/external/metacyc/metaCycRxns.mat');
% load KEGG rxns to get rxnName
load('/external/kegg/keggRxns.mat')
cd(current_path)

% load Modelseed rxn and change rev according to that
try
    websave('reactions.tsv','https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/dev/Biochemistry/reactions.tsv');
catch
    warning('reactions.tsv was not successfully downloaded, check if directory for reactions.tsv has changed on github.com/ModelSEED/ModelSEEDDatabase/Biochemistry');
end
fid2 = fopen('reactions.tsv');
format = repmat('%s ',1,23);
format = strtrim(format);
rxn_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
fclose(fid2);
for i = 1:length(rxn_temp)
    SEEDrxn(:,i) = rxn_temp{i};
end

% Map to seed part
rxn_query = intersect(MNXrxn(:,1),MNXrefRxns.SEEDMNXid); % find intersection of MNXrxn and mappingIDof mnx and metaNet
[~,rxnIdx] = ismember(rxn_query,MNXrefRxns.SEEDMNXid); % index those rxns in ref
rxn_query = MNXrefRxns.SEEDxref(rxnIdx); % find coresponding metacyc IDs
rxn_query = intersect(rxn_query,SEEDrxn(:,1)); % intersect with metacyc rxns with equation and rev
[~,rxnIdx] = ismember(rxn_query,SEEDrxn(:,1)); % index those rxns
rxn_query = SEEDrxn(rxnIdx(~strcmp(SEEDrxn(rxnIdx,9),'=')& ~strcmp(SEEDrxn(rxnIdx,9),'?')),1); %  CHOOSE REVERSIVERBILITY IN SEED rxn database
seed_equ = SEEDrxn(rxnIdx(~strcmp(SEEDrxn(rxnIdx,9),'=')& ~strcmp(SEEDrxn(rxnIdx,9),'?')),7);
seed_rev = SEEDrxn(rxnIdx(~strcmp(SEEDrxn(rxnIdx,9),'=')& ~strcmp(SEEDrxn(rxnIdx,9),'?')),9);
seed_rev = strrep(seed_rev,'<','-1');
seed_rev = strrep(seed_rev,'>','1');
seed_equ = strrep(seed_equ,' <= (',' <=> (');
seed_equ = strrep(seed_equ,'(','');
seed_equ = strrep(seed_equ,')','');
[~,rxnIdx] = ismember(rxn_query,MNXrefRxns.SEEDxref);
rxn_query = MNXrefRxns.SEEDMNXid(rxnIdx);
[~,rxnIdx] = ismember(rxn_query,MNXrxn(:,1));
MNX_equ = MNXrxn(rxnIdx,2);

% Construct S
[S_MNX, mets_MNX]=constructS(MNX_equ);
[S_seed, mets_seed]=constructS(seed_equ);
mets_MNX = strrep(mets_MNX,'@MNXD2','');
mets_MNX = strrep(mets_MNX,'@MNXD1','');
mets_seed = strrep(mets_seed,'[0]','');
mets_seed = strrep(mets_seed,'[1]','');
% map mets_metacyc into MNX IDs
mets_mapped_seed = mapIDsViaMNXref('mets',mets_seed,'SEED','MetaNetX');
for i = 1:length(rxn_query)
    mets_seed_temp = find(S_seed(:,i));
    mets_seed_coef = S_seed(mets_seed_temp,i);
    mets_seed_temp = mets_mapped_seed(mets_seed_temp);
    mets_seed_coef = mets_seed_coef(cellfun(@(x)~isempty(x(:)),mets_seed_temp));
    mets_seed_temp = mets_seed_temp(cellfun(@(x)~isempty(x(:)),mets_seed_temp));
    mets_MNX_temp = find(S_MNX(:,i));
    mets_MNX_coef = S_MNX(mets_MNX_temp,i);
    mets_MNX_temp = mets_MNX(mets_MNX_temp);
    [~,idx] = ismember(mets_seed_temp,mets_MNX_temp);
    mets_seed_coef = mets_seed_coef(find(idx)); % find coef which idx is not empty
    prod = mets_MNX_coef(idx(idx~=0)).* mets_seed_coef;
    if all(prod(:)<0)
        rev(i,1) = num2cell(-1*str2double(seed_rev(i,1))); % stand for reverse the reactants with the product
    else
        rev(i,1) = num2cell(str2double(seed_rev(i,1)));
    end
end
% update the rev information
MNXrxn(:,7) = {0};
MNXrxn(:,8) = {''};
[~,rxnIdx] = ismember(rxn_query,MNXrxn(:,1)); % index those rxns
MNXrxn(rxnIdx,7) = rev; % 0 stand for <--> 1 stand for -> -1 stand for <-
MNXrxn(rxnIdx,8) = {'rev from SEED'};

% Index the rxnID with MNXrxn
% Map to metacyc part
rxn_query = intersect(MNXrxn(:,1),MNXrefRxns.MetaCycMNXid); % find intersection of MNXrxn and mappingIDof mnx and metaNet
[~,rxnIdx] = ismember(rxn_query,MNXrefRxns.MetaCycMNXid); % index those rxns in ref
rxn_query = MNXrefRxns.MetaCycxref(rxnIdx); % find coresponding metacyc IDs
rxn_query = intersect(rxn_query,metacycRxn.metaCycRxns.rxns); % intersect with metacyc rxns with equation and rev
[~,rxnIdx] = ismember(rxn_query,metacycRxn.metaCycRxns.rxns); % index those rxns
rxn_query = metacycRxn.metaCycRxns.rxns(rxnIdx(metacycRxn.metaCycRxns.rev(rxnIdx) == 0)); %
seed_equ = metacycRxn.metaCycRxns.equations(rxnIdx(metacycRxn.metaCycRxns.rev(rxnIdx) == 0));
[~,rxnIdx] = ismember(rxn_query,MNXrefRxns.MetaCycxref);
rxn_query = MNXrefRxns.MetaCycMNXid(rxnIdx);
[~,rxnIdx] = ismember(rxn_query,MNXrxn(:,1));
MNX_equ = MNXrxn(rxnIdx,2);

% Construct S
[S_MNX, mets_MNX]=constructS(MNX_equ);
[S_metacyc, mets_metacyc]=constructS(seed_equ);
mets_MNX = strrep(mets_MNX,'@MNXD2','');
mets_MNX = strrep(mets_MNX,'@MNXD1','');
% map mets_metacyc into MNX IDs
mets_mapped_metacyc = mapIDsViaMNXref('mets',mets_metacyc,'MetaCyc','MetaNetX');
for i = 1:length(rxn_query)
    mets_seed_temp = find(S_metacyc(:,i));
    mets_seed_coef = S_metacyc(mets_seed_temp,i);
    mets_seed_temp = mets_mapped_metacyc(mets_seed_temp);
    mets_seed_coef = mets_seed_coef(cellfun(@(x)~isempty(x(:)),mets_seed_temp));
    mets_seed_temp = mets_seed_temp(cellfun(@(x)~isempty(x(:)),mets_seed_temp));
    mets_MNX_temp = find(S_MNX(:,i));
    mets_MNX_coef = S_MNX(mets_MNX_temp,i);
    mets_MNX_temp = mets_MNX(mets_MNX_temp);
    [~,idx] = ismember(mets_seed_temp,mets_MNX_temp);
    mets_seed_coef = mets_seed_coef(find(idx)); % find coef which idx is not empty
    prod = mets_MNX_coef(idx(idx~=0)).* mets_seed_coef;
    if all(prod(:)<0)
        rev(i,1) = num2cell(-1); % stand for reverse the reactants with the product
    else
        rev(i,1) = num2cell(1);
    end
end
% update the rev information
[~,rxnIdx] = ismember(rxn_query,MNXrxn(:,1)); % index those rxns
MNXrxn(rxnIdx,7) = rev; % 0 stand for <--> 1 stand for -> -1 stand for <-
MNXrxn(rxnIdx,8) = {'rev from Metacyc'};

% find kegg rxn id
rxns_mapped_kegg = mapIDsViaMNXref('rxns',MNXrxn(:,1),'MetaNetX','KEGG');
MNXrxn(:,9) = rxns_mapped_kegg;
%% add rxnName for all rxns
%seed part
rxns_mapped_SEED = mapIDsViaMNXref('rxns',MNXrxn(:,1),'MetaNetX','SEED');
rxn_SEED = MNXrxn(cellfun(@(x)~isempty(x(:)),rxns_mapped_SEED),1);
rxns_mapped_SEED = rxns_mapped_SEED(cellfun(@(x)~isempty(x(:)),rxns_mapped_SEED),1);
[~,rxnIdx] = ismember(rxns_mapped_SEED,SEEDrxn(:,1)); % index those rxns
rxn_SEED_name = SEEDrxn(rxnIdx,3);
[~,rxnIdx] = ismember(rxn_SEED,MNXrxn(:,1)); % index those rxns
MNXrxn(:,10) = MNXrxn(:,1);
MNXrxn(rxnIdx,10) = rxn_SEED_name;

% kegg part
rxns_mapped_kegg = mapIDsViaMNXref('rxns',MNXrxn(:,1),'MetaNetX','KEGG');
rxn_KEGG = MNXrxn(cellfun(@(x)~isempty(x(:)),rxns_mapped_kegg),1);
rxns_mapped_kegg = rxns_mapped_kegg(cellfun(@(x)~isempty(x(:)),rxns_mapped_kegg),1);
[~,rxnIdx] = ismember(rxns_mapped_kegg,model.rxns); % index those rxns load kegg model? model.rxnsstands for all rxnIDs for KEGG rxns
rxn_KEGG_name = model.rxnNames(rxnIdx(rxnIdx~=0));
rxn_KEGG = rxn_KEGG(find(rxnIdx~=0));
[~,rxnIdx] = ismember(rxn_KEGG,MNXrxn(:,1)); % index those rxns
MNXrxn(rxnIdx,10) = rxn_KEGG_name;
save('../../data/databases/MNXrxn.mat','MNXrxn')

%% met part

fid2 = fopen('/Users/feiranl/Downloads/chem_prop.tsv'); %MNX reactions databse with equation
format = repmat('%s ',1,9);
format = strtrim(format);
met_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(met_temp)
    MNXmet(:,i) = met_temp{i};
end
commentLines = startsWith(MNXmet(:,1),'#');
MNXmet(commentLines,:) = []; % MNX_ID	Description	Formula	Charge	Mass	InChI	SMILES	Source	InChIKey
fclose(fid2);

% find kegg id
mets_mapped_kegg = mapIDsViaMNXref('mets',MNXmet(:,1),'MetaNetX','KEGG');
MNXmet(:,10) = mets_mapped_kegg;

% find chebi id
mets_mapped_chebi = mapIDsViaMNXref('mets',MNXmet(:,1),'MetaNetX','ChEBI');
MNXmet_temp = MNXmet(cellfun(@(x)~isempty(x(:)),mets_mapped_chebi),1);
mets_mapped_chebi = mets_mapped_chebi(cellfun(@(x)~isempty(x(:)),mets_mapped_chebi),1);
aaa = repmat({'chebi'},length(mets_mapped_chebi),1);
aaa = join([aaa,mets_mapped_chebi],':');
[~,rxnIdx] = ismember(MNXmet_temp,MNXmet(:,1)); % index those rxns
MNXmet(:,11) = {''};
MNXmet(rxnIdx,11) = cellstr(aaa);
rxnIdx = startsWith(MNXmet(:,8),'chebi');
MNXmet(rxnIdx,11) = MNXmet(rxnIdx,8);
save('../../data/databases/MNXmet.mat','MNXmet')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function downloadMNXdb(files,mnxPath)
%downloadMNXdb  Download MetaNetX database files if they cannot be found.
%
%   files       string or cell array of strings, with the MetaNetX files to
%               be downloaded. Options: 'chem_xref', 'chem_prop',
%               'reac_xref', 'reac_prop' or 'all'.
%   mnxPath     string of path where MetaNetX reference files are
%               stored. (opt, default to RAVENdir/external/metanetx) To
%               download to current folder, specify pwd().
%
% Usage: downloadMNXdb(files,mnxPath)
%
% Eduard Kerkhoven, 2020-05-04

if nargin<2
    [ST, I]=dbstack('-completenames');
    mnxPath=fileparts(fileparts(fileparts(ST(I).file)));
    mnxPath=fullfile(mnxPath,'external','metanetx');
end

if ischar(files)
    if strcmp(files,'all')
        files={'chem_xref','chem_prop','reac_prop','reac_xref'};
    else
        files={files};
    end
end

for k=1:length(files)
    if ~exist(fullfile(mnxPath,[files{k},'.tsv']), 'file')
        fprintf('File %s.tsv cannot be found and will be downloaded from MetaNetX.org.\n',files{i});
        websave(fullfile(mnxPath,[files{k},'.tsv']),...
            ['https://www.metanetx.org/cgi-bin/mnxget/mnxref/',files{k},'.tsv']);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slashPos = getSlashPos(path)
slashPos = strfind(path,'\');       %Windows
if isempty(slashPos)
    slashPos = strfind(path,'/');   %MAC/Linux
end
end

