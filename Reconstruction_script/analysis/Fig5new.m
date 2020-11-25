%% Substrate evolution analysis
%% load data
current_path = pwd;

load('../ModelFiles/model_original_new.mat')
% based on clade
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
Strain_information(:,i) = temp{i};
end
fclose(fid2);
clades = unique(Strain_information(:,2));
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
clades(1) = [];

fid2 = fopen('../data/physiology/Biolog_substrate.tsv');
format = repmat('%s ',1,333);
format = strtrim(format);
substrate = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(substrate)
    data(:,i) = substrate{i};
end
strainlist = data(1,5:end);
data(1,:) = [];
SubBiologName = data(:,1);
SubModelName = data(:,2);
data(:,1:4) = [];
fclose(fid2);

% replace the experimental data with numbers
data = strrep(data,'n','nan');
data = strrep(data,'v','1');
data = cellfun(@str2num, data, 'UniformOutput', false);

load('../Reconstruction/otherchanges/FBAresultnew2.mat')
FBAresult(abs(FBAresult)>0) = 1;
data(isnan(cell2mat(data))) = num2cell(FBAresult(isnan(cell2mat(data))));
clear fid2 Subtype SubCondition format temp substrate FBAresult;

%% special/general enzyme for each rxn
metname_nocomp = model_original.metNames;
for i=1:numel(model_original.comps)
    metname_nocomp=regexprep(metname_nocomp,['\[', model_original.compNames{i}, '\]$'],'');
end
[metname_nocomp_unique,~,b] = unique(metname_nocomp,'stable');
% remake the S matrix
newmatrix = zeros(length(metname_nocomp_unique),length(model_original.rxns));
for i = 1:length(model_original.rxns)
    idx = find(model_original.S(:,i));
    newmatrix(b(idx),i) = model_original.S(idx,i);
end
% find the rxn with grRules
idx = ismember(model_original.grRules,'');
newmatrix = newmatrix(:,~idx);
rxngenemat = model_original.rxnGeneMat(~idx,:);
[newmatrix_unique,a,b] = unique(newmatrix','rows');
for i = 1:length(newmatrix_unique)% find all grRules related to the same reaction in diff comps
    idx = find(b == i);
    if length(idx) > 1
        newrxngenemat(i,:) = sum(rxngenemat(idx,:),1);%union the genes
    else
        newrxngenemat(i,:) = rxngenemat(idx,:);
    end
end
% get the count for each OG
OG_rxn_count = [model_original.genes,num2cell(sum(newrxngenemat,1))'];

cdfplot(cell2mat(OG_rxn_count(:,2)))
ylabel('Percentage','FontSize',10,'FontName','Helvetica','Color','k');
xlabel('OG linked rxn number','FontSize',10,'FontName','Helvetica','Color','k');

%[~,~,index] = xlsread('../data/substrateUsageGene.xlsx','index');
[~,~,index] = xlsread('../data/substrateUsageGene_35.xlsx','index');
rxn = index(startsWith(index(:,1),'r_'),1);
rxn_tmp = join(rxn,';');
rxn_tmp = split(rxn_tmp,';');
[~,idx] = ismember(rxn_tmp,model_original.rxns);
OG_all = model_original.genes(find(sum(model_original.rxnGeneMat(idx,:),1)));
[~,idx] = ismember(OG_all,OG_rxn_count(:,1));
OG_all_general = OG_all(find(cell2mat(OG_rxn_count(idx,2))>1));
clear i newrxngenemat idx newmatrix_unique b newmatrix metname_nocomp_unique metname_nocomp idx

%% Fig2a substrate usage in each clade/ metabolic innovation /loss
k = 0;
trait_total = {};
trait = {};
step = [];

for i = 1:length(index(:,1))
    i
    if endsWith(index{i,1},'_model')
        trait_tmp = replace(index{i,1},'_model','');
         k = 0;
    elseif endsWith(index{i,1},'_exp')
        k = 0;
        trait_tmp = replace(index(i,1),'_exp','');
        trait_total = [trait_total;trait_tmp];
    else
        if startsWith(index{i,1},'r_')
            if isnan(index{i,6}) 
                k = k + 1;
            elseif startsWith(index{i,6},'route')
                k = 1;
            end
        
        step = [step;k];
        trait = [trait;trait_tmp];
        end
    end
end
index = index(startsWith(index(:,1),'r_'),:);

%get strains_sort_clade ready
group = [];
strains_sortclade = [];
for i = 1:length(clades)
idx = ismember(Strain_information(:,2),clades(i));
[~,ID] = ismember(Strain_information(idx,1),strainlist);
group = [group;i*ones(length(ID(ID~=0)),1)];
strains_sortclade = [strains_sortclade;strainlist(ID(ID~=0))'];
end

clear trait_tmp temp idx ID genes gene_counts x i data_selected count fid2 k a ans

%% define metabolic innovation/loss
% BYCA definition
fid2 = fopen('../data/physiology/BYCA_phenotype.tsv');
format = repmat('%s ',1,6);
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
    BYCA(:,i) = temp{i};
end

% sort BYCA based on trait
[~,BYCAidx] = ismember(trait_total,BYCA(:,2));
BYCA = BYCA(BYCAidx,:); % sort as trait_total

[~,BYCAidx] = ismember(trait_total,SubModelName);
data = data(BYCAidx,:); %sort as trait_total

for i = 1:length(strains_sortclade)
    [~,ID] = ismember(strains_sortclade(i,1),strainlist);
    data_selected = data(:,ID);
    clade_av(i) = sum(cell2mat(data_selected),1,'omitnan');
    result = cell2mat(data_selected) - cell2mat(cellfun(@str2num, BYCA(BYCAidx~=0,6), 'UniformOutput', false));
    gain(i) = length(find(result > 0.85));
    gain_result(i) = join(BYCA(result > 0.85,1),';');
    loss(i) = length(find(result < -0.85));
    loss_result(i) = join(BYCA(result < -0.85,1),';');
end
%% Gene expansion/contraction for each clade
fid2 = fopen('../data/gene_expansion_contraction_subtrate.tsv');
format = '%s %s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
all_gene_change(:,i) = temp{i};
end
fclose(fid2);
[~,idx] = ismember(lower(strains_sortclade),lower(all_gene_change(:,1)));
gene_expansion = cellfun(@str2num, all_gene_change(idx,2), 'UniformOutput', false);
gene_contraction = cellfun(@str2num, all_gene_change(idx,3), 'UniformOutput', false);

fid = fopen('../data/aaa.json', 'r');
data = textscan(fid,'%s','TreatAsEmpty','NA');
data = join(data{1,1}','');
exp_con = jsondecode(data{1});
fid = fclose(fid);

format = '%s %s';
fID       = fopen('../../find_homolog_for_panID_py/data/representatives.tsv');
mapping  = textscan(fID,format,'Delimiter','\t','HeaderLines',1);
mappingID.OGIDs      = mapping{1};
mappingID.panIDs = mapping{2};
fclose(fID);

% add back the OG ID to the panmodel
[~,ID] = ismember(model_original.genes,strrep(mappingID.panIDs,'Saccharomyces_cerevisiae@',''));
model_original.OGIDs = mappingID.OGIDs(ID);
mappingID.panIDs = strrep(mappingID.panIDs,'Saccharomyces_cerevisiae@','');
exp_sub = zeros(length(rxn),length(strains_sortclade));
con_sub = zeros(length(rxn),length(strains_sortclade));
exp_sub_count = zeros(length(strains_sortclade),1);
con_sub_count = zeros(length(strains_sortclade),1);
for i = 1:length(strains_sortclade)
    
    disp(['No. ',num2str(i), ': ',strains_sortclade{i}])
    [~,idx] = ismember(strains_sortclade(i),{exp_con.organism});
    tmp = exp_con(idx).expansion;
    tmp2 = exp_con(idx).contraction;
    cd(inputpath)
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    if ~isempty(tmp)
        [~,idx2] = ismember(tmp,mappingID.OGIDs);
        tmp = mappingID.panIDs(idx2); % match model.proteins
        [~,idx3] = ismember(tmp,model.proteins);
        exp_sub_count(i) = length(idx3(idx3~=0));
        exp_rxnlst = model.rxns(find(sum(model.rxnGeneMat(:,idx3(idx3~=0)),2)));
        if ~isempty(exp_rxnlst)
            for j = 1:length(exp_rxnlst)
                idx4 = find(contains(rxn,exp_rxnlst(j)));
                exp_sub(idx4(idx4~=0),i)= 1;
                
            end
        end
    end
    if ~isempty(tmp2)
        [~,idx2] = ismember(tmp2,mappingID.OGIDs);
        tmp2 = mappingID.panIDs(idx2); % match model.proteins
        [~,idx3] = ismember(tmp2,model.proteins);
        con_sub_count(i) = length(idx3(idx3~=0));
        con_rxnlst = model.rxns(find(sum(model.rxnGeneMat(:,idx3(idx3~=0)),2)));
        if ~isempty(con_rxnlst)
            for j = 1:length(con_rxnlst)
                idx4 = find(contains(rxn,con_rxnlst(j)));
                con_sub(idx4(idx4~=0),i)= 1;
            end
        end
    end
end
clear model idx idx2 tmp tmp2;
%% HGT plot based on clade
cd(current_path)
fid2 = fopen('../data/substrate_bacteria_HGT.tsv');
format = '%s %s %s';
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
HGT1(:,i) = temp{i};
end
fclose(fid2);

fid2 = fopen('../data/substrate_fungi_HGT.tsv');
format = '%s %s %s';
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
HGT2(:,i) = temp{i};
end
fclose(fid2);

fid2 = fopen('../data/transporter_bacteria_HGT.tsv');
format = '%s %s %s';
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
HGT3(:,i) = temp{i};
end
fclose(fid2);

fid2 = fopen('../data/transporter_fungi_HGT.tsv');
format = '%s %s %s';
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
HGT4(:,i) = temp{i};
end
fclose(fid2);

HGT = [HGT1;HGT2;HGT3;HGT4];
HGT = split(unique(join(HGT,';'),'stable'),';');
% sort HGT by species
temp = split(HGT(:,2),'@'); % split by the gene
HGT(:,4) = temp(:,1); % species
HGT_rxn = zeros(length(rxn),length(strains_sortclade));
for i = 1:length(strains_sortclade)
    i
    HGT_rxn_temp = cell(0,1);
    ID = ismember(HGT(:,4),strains_sortclade(i,1));
    HGT_species(i) = length(HGT(ID,3));
    HGT_enzyme(i) = join(HGT(ID,2),';');
    cd(inputpath)
    if ~isempty(HGT(ID,4))
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    [~,idx] = ismember(HGT(ID,2),model.genes);
    HGT_temp = model.rxns(logical(sum(full(model.rxnGeneMat(:,idx(idx~=0))==1),2)));
    for j = 1:length(HGT_temp)
    idx_tmp = find(contains(rxn,HGT_temp(j)));
        HGT_rxn(idx_tmp,i) = 1; % which is for later use
    end
    end
end
for i = 1:length(strains_sortclade)
    i
    HGT_enzyme_temp = cell(0,1);
    ID = ismember(HGT(:,4),strains_sortclade(i,1));
    [~,idx] = ismember(HGT(ID,1),model_original.OGIDs);
    HGT_enzyme_temp = model_original.genes(idx(idx~=0));
    HGT_count(i) = length(intersect(OG_all,HGT_enzyme_temp));
end
cd(current_path)
clear temp model;

%% generalist/special list
for i = 1:length(strains_sortclade)
    disp(['No.',num2str(i),': ',strains_sortclade{i}])
       cd(inputpath)
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    for j = 1:length(rxn)
      rxn_temp = split(rxn{j},';');
        [~,idx] = ismember(rxn_temp,model.rxns);
        genes = sum(model.rxnGeneMat(idx(idx~=0),:),1);
        [~,idx] = ismember(unique(model.proteins(find(genes))),OG_rxn_count(:,1)); % OG_rxn_count stands for the rxn number mapped to each panID
        gene_counts = cell2mat(OG_rxn_count(idx(idx~=0),2)); % gene_counts stands for rxn number correponded to this gene
        if ~isempty(idx==0) % deal with specific species genes.
            genes_temp = model.genes(find(genes));
            [~,id] = ismember(genes_temp(idx == 0),model.genes);
            temp = sum(model.rxnGeneMat(:,id),1);
            gene_counts(idx == 0) = temp;
        end
        general_rxn(j,i) = sum(gene_counts > 1); % gene number with muliple reactions
    end
end

for i = 1:length(strains_sortclade)
    disp(['No.',num2str(i),': ',strains_sortclade{i}])
       cd(inputpath)
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    general_OG_count(i) = length(intersect(model.proteins,OG_all_general));
end
%% Fig1a
subplot(1,6,1)

%f1 = scatter(x_total,clade_av,15,'filled','k');f1.MarkerFaceAlpha = 0.4;
h1 = boxplot(clade_av,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
set(h1,{'linew'},{1});
camroll(-90)
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('No. traits','FontSize',12,'FontName','Helvetica');
subplot(1,6,2)
%f2 = scatter(x_total,gain,15,'filled','k');f2.MarkerFaceAlpha = 0.4;

h2 = boxplot(gain,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[255,127,0]/255);
set(h2,{'linew'},{1});
camroll(-90)
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('gained','FontSize',12,'FontName','Helvetica');
subplot(1,6,3)

%f3 = scatter(x_total,loss,15,'filled','k');f3.MarkerFaceAlpha = 0.4;
h3 = boxplot(loss,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[49,163,84]/255);
set(h3,{'linew'},{1});
camroll(-90)
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('lost ','FontSize',12,'FontName','Helvetica');
subplot(1,6,4)

%f4 = scatter(x_total,HGT_species,15,'filled','k');f4.MarkerFaceAlpha = 0.4;
h4 = boxplot(HGT_count,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[215,48,39]/255);
set(h4,{'linew'},{1});
camroll(-90)
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Substrate HGT','FontSize',12,'FontName','Helvetica');
subplot(1,6,5)

%f5 = scatter(x_total,gene_expansion,15,'filled','k');f5.MarkerFaceAlpha = 0.4;
h5 = boxplot(exp_sub_count,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[118,42,131]/255);
set(h5,{'linew'},{1});
camroll(-90)
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Gene expansion','FontSize',12,'FontName','Helvetica');
%subplot(1,6,6)

% %f6 = scatter(x_total,gene_contraction,15,'filled','k');f5.MarkerFaceAlpha = 0.4;
% h6 = boxplot(sum(con_sub,1),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[118,42,131]/255);
% set(h6,{'linew'},{1});
% camroll(-90)
% set(gca,'FontSize',10,'FontName','Helvetica');
% ylabel('Gene contraction','FontSize',12,'FontName','Helvetica');

subplot(1,6,6)
%f6 = scatter(x_total,sum(general_rxn,1),15,'filled','k');f6.MarkerFaceAlpha = 0.4;
h6 = boxplot(general_OG_count,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[221,28,119]/255);
set(h6,{'linew'},{1});
camroll(-90)
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('General list','FontSize',12,'FontName','Helvetica');


%% write table out for R plot
table_cor = table(clade_av',gain',loss',HGT_count',exp_sub_count,general_OG_count');
table_cor.Properties.VariableNames = {'No.traits' 'Gain' 'Loss' 'Substrate HGT','Gene expansion','General list'};
writetable(table_cor,'table.txt')
 %% Fig1d for HGT/gene_expansion/Generalist in metabolic gain
for i = 1:length(strains_sortclade)
     gain_idx = split(gain_result(1,i),';');
     gain_idx = replace(gain_idx,SubBiologName(8:end),SubModelName(8:end));
     gain_idx = intersect(gain_idx,SubModelName);
     gain_idx(strcmp(gain_idx,'')) = [];
     HGT_idx = unique(trait(find(HGT_rxn(:,i)>=1)));
     exp_idx = unique(trait(find(exp_sub(:,i)>=1)));
     GL_idx = unique(trait(general_rxn(:,i)>=1));
     contribution(i,1) = length(intersect(gain_idx,HGT_idx)); %  HGT
     contribution(i,2) = length(intersect(gain_idx,exp_idx)); %  GD
     contribution(i,3) = length(intersect(gain_idx,GL_idx)); %  GL
     contribution(i,4) = length(setdiff(gain_idx,union(union(HGT_idx,exp_idx),GL_idx))); % no reason found
end
ylables = {'HGT','GE','GL','NO REASON'};
 

% second figure
for i = 1:length(clades)
    idx = ismember(group,i);
    contribution_clade(i,1) = sum(contribution(idx,1),1);
    contribution_clade(i,2) = sum(contribution(idx,2),1);
    contribution_clade(i,3) = sum(contribution(idx,3),1);
    contribution_clade(i,4) = sum(contribution(idx,4),1);
end

% Fig 1d
h1 = bar(sum(contribution_clade,1)/sum(contribution_clade(:)),'FaceColor',[255,127,0]/255,'FaceAlpha',0.3,'EdgeColor',[255,127,0]/255,'LineWidth',1);
set(gca,'XTick',1:1:4);
set(gca,'XTickLabel',{'HGT','Gene expansion','Generalist','Others'});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Ratio of Evolution events cocurring with subtrate gain','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(0);

% for all HGT where they locate
HGT_rxn_sum = sum(HGT_rxn,2); % refer to each rxn HGT
for i = 1:max(step)
    HGT_each_step(i) = sum(HGT_rxn_sum(step==i,:)); % hGT haappens in the first translocation step; change to 2 is the first degradation rxn
end
%% Figure for HGT num based on clade
% All HGT
% fid2 = fopen('../data/all_HGT_1840.tsv');
% format = repmat('%s ',1,11);
% format = strtrim(format);
% temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
% for i = 1:length(temp)
% all_HGT(:,i) = temp{i};
% end
% fclose(fid2);
% all_HGT = all_HGT(find(strcmp(all_HGT(:,11),'Yes')),:);
% temp = split(all_HGT(:,2),'@'); % split by the gene
% all_HGT(:,3) = temp(:,1);
% temp = split(all_HGT(:,6),'-');
% all_HGT(:,12) = temp(:,4);
% origin = {'Actinoba';'Amoebozoa';'Apicomple';'BacteroC';'Chlorofl';'Cyanobacte';'Euglenozo';'Firmicut';'Gemmatimo';'Metazoa';'Proteoba';'Streptoph';'Thaumarc';'other_Ba';'other_Stra'};
% for j = 1:length(strains_sortclade)
%     j
%     ID = ismember(all_HGT(:,3),strains_sortclade(j,1));
%     HGT_species(j) = length(unique(all_HGT(ID,2)));
%     HGT_average_AI(j) = mean(cellfun(@str2num, all_HGT(ID,10)));
%     for k = 1:length(origin)
%         temp = find(startsWith(all_HGT(ID,12),origin{k}));
%         origin_HGT(j,k) = length(temp);
%     end
% end
all_HGT = HGT;
origin = unique(HGT(:,3));
for j = 1:length(strains_sortclade)
    j
    ID = ismember(all_HGT(:,4),strains_sortclade(j,1));
    HGT_species(j) = length(unique(all_HGT(ID,2)));
    for k = 1:length(origin)
        temp = find(startsWith(all_HGT(ID,3),origin{k}));
        origin_HGT(j,k) = length(temp);
    end
end

% figure
% subplot(2,1,2)
% hold on
% h2 = boxplot(HGT_species,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
% f2 = scatter(x_total,HGT_species,20,'k','filled');f2.MarkerFaceAlpha = 0.4;
% ylim([0,100])
% set(gca,'FontSize',10,'XTickLabelRotation',90)
% ylabel('All HGTs in each clade','FontSize',12,'FontName','Helvetica');
% subplot(2,1,1)
% h2 = boxplot(HGT_species,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
% f2 = scatter(x_total,HGT_species,20,'k','filled');f2.MarkerFaceAlpha = 0.4;
% ylim([200,250])
% ylabel('All HGTs in each clade','FontSize',12,'FontName','Helvetica');

%% Fig2b ALL ORIGIN OF hgt
temp = sum(origin_HGT,1);
[R C]=sort(temp,'descend');
lable = [origin(C(1:5));'Others'];

h1 = bar([R(1:5),sum(R(6:end))],'FaceColor',[56,108,176]/255,'FaceAlpha',0.3,'EdgeColor',[56,108,176]/255,'LineWidth',1);
set(gca,'XTick',1:1:6);
set(gca,'XTickLabel',lable);
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('HGT origin','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(0);

% origin_HGT_final = [origin_HGT(:,C(1:5)),sum(origin_HGT(:,C(6:end)),2)];
% for k = 1:6
% subplot(1,6,k)
% h = boxplot(origin_HGT_final(:,k),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
% set(h,{'linew'},{1});
% camroll(-90)
% ylim([0,20])
% ylabel(strrep(lable{k},'_',' '),'FontSize',12,'FontName','Helvetica');
% end

% transport HGT
% based on transport annotation from the gene annotation of eggnog and kegg
% [~,~,transporter] = xlsread('../data/transpoter_annotation_of_all_OG.xlsx');
% mappingID.OGIDs      = mapping{1};
% mappingID.panIDs = mapping{2};
% [~,ID] = ismember(HGT(:,1),mappingID.OGIDs);
% HGT(:,14) = mappingID.panIDs(ID);
% HGT_transporter = intersect(HGT(:,14),strrep(transporter(:,1),'Saccharomyces_cerevisiae@',''));
% for i = 1:length(HGT_transporter)
% num = ismember(HGT(:,14),HGT_transporter(i));
% HGT_transporter(i,2) = num2cell(sum(num));
% end
% [~,ID] = ismember(HGT_transporter(:,1),transporter(:,1));
% HGT_transporter(ID~=0,3) = transporter(ID(ID~=0),end);

% transporter HGT combine transporter and also the rxn
HGT_transporter = [HGT3;HGT4]; % annotation from gapse transporter
HGT_sub = [HGT1;HGT2]; % from modeL, bu tthey also contain transporter, so we should seperate them.
transrxn = index(strcmp(index(:,end),'trans'));
transrxn = split(join(transrxn,';'),';');
[~,idx] = ismember(transrxn,model_original.rxns);
transOG = model_original.OGIDs(find(sum(full(model_original.rxnGeneMat(idx,:)),1)));

% extracelluar enzyme
exrxn = index(strcmp(index(:,end),'extracelluar'));
exrxn = split(join(exrxn,';'),';');
[~,idx] = ismember(exrxn,model_original.rxns);
exOG = model_original.OGIDs(find(sum(full(model_original.rxnGeneMat(idx,:)),1)));


[~,idx] = ismember(HGT_sub(:,1),transOG);
HGT_transporter = [HGT_transporter;HGT_sub(idx~=0,:)];

[~,idx2] = ismember(HGT_sub(:,1),exOG);
HGT_ex = HGT_sub(idx2~=0,:);

HGT_sub = HGT_sub((idx == 0 & idx2 == 0),:);
HGT_sub = split(unique(join(HGT_sub,';'),'stable'),';');
HGT_transporter = split(unique(join(HGT_transporter,';'),'stable'),';');
HGT_ex = split(unique(join(HGT_ex,';'),'stable'),';');
h1 = bar([length(HGT_ex(:,1)),length(HGT_transporter(:,1)),length(HGT_sub(:,1))],'FaceColor',[56,108,176]/255,'FaceAlpha',0.3,'EdgeColor',[56,108,176]/255,'LineWidth',1);
set(gca,'XTick',1:1:3);
set(gca,'XTickLabel',{'Extracelluar enzymes','Transporter','cellular enzymes'});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('HGT classification','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(0);

cd(currennt_path)

%% substrate loss
% define the metabolic loss as critical enzyme loss
%result_table = SubstrateUsageGene_Table;% or load directly
[~,~,result_table] = xlsread('../data/substrateUsageGene_35.xlsx','RESULTTABLE');
result_table = result_table(10:end,2:end);
[~,~,index] = xlsread('../data/substrateUsageGene_35.xlsx','index');
[~,idx] = ismember(strains_sortclade,strainlist);
result_table = result_table(idx,:);
rxnexistence = result_table(:,startsWith(index(:,1),'r_'));
FBAresult = result_table(:,endsWith(index(:,1),'_model'));
sub_exp = result_table(:,endsWith(index(:,1),'_exp'));
sub_exp(strcmp(sub_exp,'n')) = {nan};
sub_exp(strcmp(sub_exp,'v')) = {1};
k = 0;
m = 0;
trait_total = {};
trait = {};
trait_idx = {};
for i = 1:length(index(:,1))
    i
    if endsWith(index{i,1},'_model')
        trait_tmp = replace(index{i,1},'_model','');
    elseif endsWith(index{i,1},'_exp')
        k = k+1;
        trait_tmp = replace(index(i,1),'_exp','');
        trait_total = [trait_total;trait_tmp];
    else
        m = m+1;
        step(m,1) = k;
    end
end
index = index(startsWith(index(:,1),'r_'),:);


for i = 1:length(strains_sortclade) % define the loss type
    loss_temp = find(cell2mat(sub_exp(i,:))==0); % lost sub
    for j = 1:length(loss_temp)
       rxn_idx = find(step(:,1)==loss_temp(j)); % fincd corresponding rxns
       criticalrxn_idx = cell2mat(index(step(:,1)==loss_temp(j),3)) > 0.83 & cell2mat(index(step(:,1)==loss_temp(j),4)) > 0.92;% based on the creteria of galactose
       rxn_idx = rxn_idx(criticalrxn_idx);
       if any(rxn_idx)
            temp = rxnexistence(i,rxn_idx);
            model_temp = FBAresult(i,loss_temp(j));
            if cell2mat(model_temp) == 0
                critical_loss(i,loss_temp(j)) = length(temp(cell2mat(temp)==0))/length(rxn_idx);
            else
                critical_noloss(i,loss_temp(j)) = length(temp(cell2mat(temp)==0));
            end
       end
    end
end

sum(cellfun(@any,num2cell(critical_loss)),2)
h = boxplot(ans,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(gca,'FontSize',10,'XTickLabelRotation',90)
ylabel('Substrate loss caused by critical rxn loss','FontSize',12,'FontName','Helvetica','Color','k');

% define the metabolic loss as non-critical enzyme loss
for i = 1:length(strains_sortclade) % define the loss type
    loss_temp = find(cell2mat(sub_exp(i,:))==0);
    for j = 1:length(loss_temp)
       rxn_idx = find(step(:,1)==loss_temp(j));
       criticalrxn_idx = cell2mat(index(step(:,1)==loss_temp(j),3)) > 0.83 & cell2mat(index(step(:,1)==loss_temp(j),4)) > 0.92;% based on the creteria of galactose
       nocriticalrxn_idx = ~criticalrxn_idx;% based on the creteria of galactose
       norxn_idx = rxn_idx(nocriticalrxn_idx);
       if any(norxn_idx)
            temp_no = rxnexistence(i,norxn_idx);
            model_temp = FBAresult(i,loss_temp(j));
            if cell2mat(model_temp) == 0
                no_critical_loss(i,loss_temp(j)) = length(find(cell2mat(temp_no)==0))/length(norxn_idx);
            else
                no_critical_noloss(i,loss_temp(j)) = length(find(cell2mat(temp_no)==0))/length(norxn_idx);
            end
       end
    end
end

sum(cellfun(@any,num2cell(no_critical_loss)),2)
h = boxplot(ans,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(gca,'FontSize',10,'XTickLabelRotation',90)
ylabel('Substrate loss caused by critical rxn loss','FontSize',12,'FontName','Helvetica','Color','k');

% define the metabolic loss as downstream pathway
downsteam = zeros(length(strains_sortclade),length(sub_exp(1,:)));
for i = 1:length(strains_sortclade) % define the loss type
    loss_temp = find(cell2mat(sub_exp(i,:))==0);
    for j = 1:length(loss_temp)
       rxn_idx = find(step(:,1)==loss_temp(j));
            temp = rxnexistence(i,rxn_idx);
            model_temp = FBAresult(i,loss_temp(j));
            if contains(index(rxn_idx,5),'route 2') % there are alternative pwys
                alternative = unique(index(rxn_idx,5));
                for k = 1:length(alternative)
                    idx = strcmp(index(rxn_idx,5),alternative(k));
                    if all(cell2mat(temp(idx)))
                        rxnex_tmp(k) = 1;
                    end
                end
                if cell2mat(model_temp) == 0 && any(rxnex_tmp)
                    i
                    j
                    k
                    downsteam(i,loss_temp(j)) = 1;
                end
            else
                if cell2mat(model_temp) == 0&& all(cell2mat(temp))
                    downsteam(i,loss_temp(j)) = 1;
                end
            end       
    end
end

sum(cellfun(@any,num2cell(downsteam)),2)
 h = boxplot(ans,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(gca,'FontSize',10,'XTickLabelRotation',90)
ylabel('Substrate loss caused by downstream pathway','FontSize',12,'FontName','Helvetica','Color','k');

% define the metabolic loss as generalist

for i = 1:length(strains_sortclade) % define the loss type
    loss_temp = find(cell2mat(sub_exp(i,:))==0);
    for j = 1:length(loss_temp)
       rxn_idx = find(step(:,1)==loss_temp(j));
            temp = rxnexistence(i,rxn_idx);
            model_temp = FBAresult(i,loss_temp(j));
            general_Temp = full(general_rxn(rxn_idx,i));
            model_temp = FBAresult(i,loss_temp(j));
            %if cell2mat(model_temp) ~= 0 && any(general_Temp) %%%% not sure whether we should define this one as model yes
            if any(general_Temp)
                generalist(i,loss_temp(j)) = length(find(general_Temp> 0))/length(rxn_idx);
            end
    end
end

sum(cellfun(@any,num2cell(generalist)),2)
 h = boxplot(ans,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(gca,'FontSize',10,'XTickLabelRotation',90)
ylabel('general list','FontSize',12,'FontName','Helvetica','Color','k');

% Fig 2e generate a figure of Loss
 for i = 1:length(strains_sortclade)
     loss_idx = split(loss_result(1,i),';');
     loss_idx = replace(loss_idx,SubBiologName(8:end),SubModelName(8:end));
     loss_idx = intersect(loss_idx,SubModelName);
     loss_idx(strcmp(loss_idx,'')) = [];
     critical_idx = unique(trait_total(find(critical_loss(i,:)>0)));
     no_critical_idx = unique(trait_total(find(no_critical_loss(i,:)> 0 )));
     dp_idx = unique(trait_total(downsteam(i,:)>0));
     general_idx = unique(trait_total(generalist(i,:)>0));
     contribution(i,1) = length(intersect(loss_idx,critical_idx)); %  critical
     contribution(i,2) = length(intersect(loss_idx,no_critical_idx)); %  no_critical
     contribution(i,3) = length(intersect(loss_idx,dp_idx)); %  downstream pathway
     contribution(i,4) = length(intersect(loss_idx,general_idx)); %  general_list
     contribution(i,5) = length(setdiff(loss_idx,union(union(critical_idx,no_critical_idx),union(dp_idx,general_idx)))); % no reason found
 end

 % contribution(:,4) = []; remove generalist condition
h1 = bar(sum(contribution,1)/sum(contribution(:)),'FaceColor',[49,163,84]/255,'FaceAlpha',0.3,'EdgeColor',[49,163,84]/255,'LineWidth',1);
set(gca,'XTick',1:1:4);
set(gca,'XTickLabel',{'Highly-correlated rxn','Non-highly correlated rxn','Downstream pathway','Others'});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('No. subtrate loss','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(0);
set(gca,'FontSize',10,'FontName','Helvetica');

clearvars general_idx dp_idx loss_idx no_critical_idx critical_idx temp rxnexistence result norxn_idx nocritical_noloss no_critical_noloss no_critical_loss loss_temp j k
clearvars reducedModel nocriticalrxn_idx m idx ID id i h h1 genes_temp genes sub_exp temp_no trait trait_1 trait_tmp trait_total format fid2 existence
%% Fig 2c critical enzyme
critical_idx = find(cell2mat(index(:,3)) > 0.83 & cell2mat(index(:,4)) > 0.92); % based on the creteria of galactose
% get trait
[~,~,index] = xlsread('../data/substrateUsageGene_35.xlsx','index');
rxn = index(startsWith(index(:,1),'r_'),1);
for i = 1:length(index(:,1))
    if endsWith(index{i,1},'_exp')
        trait_1 = replace(index{i,1},'_exp','');
        trait{i} = trait_1;
    else
        trait{i} = trait_1;
    end
end
trait_total = unique(trait);
h = bar([length(unique(trait(critical_idx))),length(trait_total)-length(unique(trait(critical_idx)))],'FaceColor','k','FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
set(gca,'FontSize',10,'FontName','Helvetica');
%ylim([0 20])
xticks([1,2]);
xticklabels({'Highly-correlated rxn','Non-highly correlated rxn'})
ylabel('No. substrates','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(0);

%% Fig2f false positives
[~,~,index] = xlsread('../data/substrateUsageGene_35.xlsx','index');
[~,idx] = ismember(strains_sortclade,strainlist);
result_table = result_table(idx,:); % sort the result table based on the strains_sortclade
rxnexistence = result_table(:,startsWith(index(:,1),'r_'));
FBAresult = result_table(:,endsWith(index(:,1),'_model'));
sub_exp = result_table(:,endsWith(index(:,1),'_exp'));
sub_exp(strcmp(sub_exp,'n')) = {nan};
sub_exp(strcmp(sub_exp,'v')) = {1};
k = 0;
m = 0;
trait_total = [];
for i = 1:length(index(:,1))
    i
    if endsWith(index{i,1},'_model')
        trait_tmp = replace(index{i,1},'_model','');
    elseif endsWith(index{i,1},'_exp')
        k = k+1;
        trait_tmp = replace(index(i,1),'_exp','');
        trait_total = [trait_total;trait_tmp];
    else
        m = m+1;
        step(m,1) = k;
    end
end

index = index(startsWith(index(:,1),'r_'),:);

for i = 1:length(strains_sortclade) % define the loss type
    nomatch = find(cell2mat(sub_exp(i,:))==0 & cell2mat(FBAresult(i,:))~=0);
    for j = 1:length(nomatch)
        rxn_idx = find(step(:,1)==nomatch(j));
        model_temp = FBAresult(i,nomatch(j));
        general_Temp = full(general_rxn(rxn_idx,i));
        if cell2mat(model_temp) ~= 0 && any(general_Temp)
            nomatch_generalist(i,nomatch(j)) = length(find(general_Temp> 0))/length(rxn_idx);
        else
            nomatch_nogeneralist(i,nomatch(j)) = 1;
        end
    end
end

for i = 1:length(strains_sortclade) % define the loss type
    nomatch = find(cell2mat(sub_exp(i,:))~=0 & cell2mat(FBAresult(i,:))==0);
    for j = 1:length(nomatch)
            nomatch_missingannotation(i,nomatch(j)) = 1;
    end
end

count_nomatch_missingannotation = sum(cellfun(@any,num2cell(nomatch_generalist)),2);
count_nomatch_generalist = sum(cellfun(@any,num2cell(nomatch_generalist)),2);
count_nomatch_nogeneralist = sum(cellfun(@any,num2cell(nomatch_nogeneralist)),2);
h = bar([sum(count_nomatch_generalist),sum(count_nomatch_nogeneralist)]/(sum(count_nomatch_generalist) + sum(count_nomatch_nogeneralist)),'FaceColor','k','FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
ylabel('No. false positive','FontSize',10,'FontName','Helvetica','Color','k');
set(gca,'FontSize',10,'XTickLabelRotation',90)
xticklabels({'Generalist','Others'})
