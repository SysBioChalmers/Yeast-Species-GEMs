%% Substrate evolution analysis
%% load data
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
group = [];
clade_av = [];

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

clear fid2 Subtype SubCondition format temp substrate

%% first figure of special/general enzyme for each substrate
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
clear i newrxngenemat idx newmatrix_unique b newmatrix metname_nocomp_unique metname_nocomp

% find all genes with substrate usage
[~,~,index] = xlsread('../data/substrateUsageGene1.xlsx','index');
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
        rxn_temp = split(index{i,1},';');
        [~,idx] = ismember(rxn_temp,model_original.rxns);
        genes = sum(model_original.rxnGeneMat(idx,:),1);
        [~,idx] = ismember(model_original.genes(find(genes)),OG_rxn_count(:,1));

            gene_counts = OG_rxn_count(idx,2);
            ratio(m,1) = length(find(cell2mat(gene_counts)<=1))/length(gene_counts);
            ratio(m,2) = k;

    end
end

ratio(:,3) = 1; % refer to one rxn
ratio(isnan(ratio(:,1)),1) = -1; % replace the empty with -1
hold on
color = redblue(101);
for i = 1:length(unique(ratio(:,2)))
    mmm = ratio(:,2) == i;
    h2 = bar(i,ratio(mmm,3),'stacked','EdgeColor','k','Facecolor','r','LineWidth',0.5);
    x = ratio(mmm,1);
    for j = 1:length(ratio(mmm,3))
        %h2(j).FaceAlpha = 0.4;
        if x(j) == -1
            h2(j).FaceColor = [189,189,189]/255;
        else
            h2(j).FaceColor = color(round(x(j)*100)+1,:);
        end
    end
end
xticks([1:1:32]);
set(gca,'XTickLabel',trait_total);
xtickangle(90);
ylabel('No. Enzymes','FontSize',20,'FontName','Helvetica','Color','k');

clear metname_nocomp h2 trait_total ratio mmm color x k m i trait_temp rxn_temp idx j 

%% FigS10 clade trait existence based on each clade and each trait- supplementary figure
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
clades(1) = [];
for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    [~,ID] = ismember(Strain_information(idx,1),strainlist);
    data_slected = data(~strcmp(SubModelName ,''),ID(ID~=0));
    for j = 1:length(data_slected(:,1))
        m(i,j)   = sum(sum(cell2mat(data_slected(j,:))==1))/length(data_slected(1,:)); % with this phenotype
        n(i,j)   = sum(sum(cell2mat(data_slected(j,:))==0))/length(data_slected(1,:)); % without phenotype
    end
end

set(gcf,'position',[200 0 1550 1350]);
lable = SubBiologName(~strcmp(SubModelName ,''));
lable = strrep(lable,' Fermentation','*');
lable = strrep(lable,'N-Acetyl-D-glucosamine','GlcNAc');
lable = strrep(lable,'"','');
for i = 1:length(data_slected(:,1))
    subplot(8,6,i);
    bar([m(:,i),n(:,i)],'EdgeColor','k','LineWidth',0.5)
    ylabel(lable{i},'FontSize',12,'FontName','Helvetica');
end

clear  lable data_selected m n j i idx ID

%% Fig2a substrate usage in each clade/ metabolic innovation /loss
% BYCA definition
fid2 = fopen('../data/physiology/BYCA_phenotype.tsv');
format = repmat('%s ',1,6);
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
    BYCA(:,i) = temp{i};
end

clade_av = [];
group = [];
x_total = [];
strains_sortclade = [];
[~,BYCAidx] = ismember(BYCA(:,1),SubBiologName);
for i = 1:length(clades)
idx = ismember(Strain_information(:,2),clades(i));
[~,ID] = ismember(Strain_information(idx,1),strainlist);
data_slected = data(BYCAidx(BYCAidx~=0),ID(ID~=0));
count = sum(cell2mat(data_slected),1,'omitnan');
clade_av = [clade_av;count'];
group = [group;i*ones(length(ID(ID~=0)),1)];
x=ones(length(ID(ID~=0)),1).*(1+(rand(length(ID(ID~=0)))-0.5)/15);
x = x(:,1).*i;
x_total = [x_total;x];
strains_sortclade = [strains_sortclade;strainlist(ID(ID~=0))'];
end
clear trait_tmp temp idx ID genes gene_counts x i data_selected count fid2


%% define metabolic innovation/loss
for i = 1:length(strains_sortclade)
    [~,ID] = ismember(strains_sortclade(i,1),strainlist);
    data_selected = data(BYCAidx(BYCAidx~=0),ID(ID~=0));
    result = cell2mat(data_selected) - cell2mat(cellfun(@str2num, BYCA(BYCAidx~=0,6), 'UniformOutput', false));
    gain(i) = length(find(result > 0.75));
    gain_result(i) = join(BYCA(result > 0.75,1),';');
    loss(i) = length(find(result < -0.75));
    loss_result(i) = join(BYCA(result < -0.75,1),';');
end

%% HGT plot based on clade
fid2 = fopen('../data/substrate_usage_new_HGT.tsv');
format = repmat('%s ',1,12);
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
HGT(:,i) = temp{i};
end
fclose(fid2);
HGT = HGT(find(strcmp(HGT(:,12),'Yes')),:);
[~,~,index] = xlsread('../data/substrateUsageGene1.xlsx','index');
rxn = index(startsWith(index(:,1),'r_'),1);
padj = index(startsWith(index(:,1),'r_'),2);
trait = index(endsWith(index(:,1),'_exp'),1);
trait = replace(trait,'_exp','');
inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';
current_path = pwd;
% sort HGT by species
temp = split(HGT(:,3),'@'); % split by the gene
HGT(:,3) = temp(:,1);
HGT_rxn = zeros(length(rxn),length(strains_sortclade));
for i = 1:length(strains_sortclade)
    i
    HGT_rxn_temp = cell(0,1);
    ID = ismember(HGT(:,3),strains_sortclade(i,1));
    HGT_species(i) = length(HGT(ID,3));
    HGT_enzyme(i) = join(HGT(ID,2),';');
    cd(inputpath)
    if length(HGT(ID,3)) ~=0
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    [~,idx] = ismember(HGT(ID,2),model.genes);
    HGT_temp = model.rxns(logical(sum(full(model.rxnGeneMat(:,idx)==1),2)));
    for j = 1:length(HGT_temp)
    idx_tmp = find(contains(rxn,HGT_temp(j)));
        HGT_rxn(idx_tmp,i) = 1; % which is for later use
    end
    end
end
cd(current_path)

%% Gene duplication for each enzyme
[~,~,index] = xlsread('../data/substrateUsageGene1.xlsx','index');
rxn = index(startsWith(index(:,1),'r_'),1);
for i = 1:length(index(:,1))
    if endsWith(index{i,1},'_exp')
        trait_1 = replace(index{i,1},'_exp','');
        trait{i} = trait_1;
    else
        trait{i} = trait_1;
    end
end
trait = replace(trait,'_exp','');
trait = trait(startsWith(index(:,1),'r_'),1);
current_path = pwd;
for i = 1:length(strains_sortclade)
    i
    cd(current_path)
    fileName    = ['../../../Multi_scale_evolution/pan_genome/result/id_mapping/',strains_sortclade{i},'.tsv'];
    fID         = fopen(fileName);
    protData    = textscan(fID,'%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    fclose(fID);
    geneID_core = protData{2}; % geneID in 343 yeast species with @seq
    panID_final = protData{5}; % panID
    panID_final = strrep(panID_final,'Saccharomyces_cerevisiae@','');
    cd(inputpath)
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    for j = 1:length(rxn)
        rxn_s = split(rxn{j},';');
        [~,ID] = ismember(rxn_s,model.rxns);
        if sum(model.rxnGeneMat(ID(ID~=0),:)) ~=0
            genelist = model.genes(logical(sum(model.rxnGeneMat(ID(ID~=0),:),1)));
            old_grRule(j,i) = join(genelist,';');
            old_grRule_count(j,i) = length(genelist);
            % index the grRule by panID_final
            [~,idx] = ismember(genelist,panID_final);
            if ~all(idx) || isempty(idx)
                warning([genelist{idx==0},'cannot be found in panID table in ', strains_sortclade{i}])
            else
                a = find(ismember(panID_final,genelist));
                new_grRule(j,i) = join(geneID_core(a),';');
                new_grRule_count(j,i) = length(panID_final(a));
            end
        end
    end
end
cd(current_path)
% 
GD_rxn = new_grRule_count - old_grRule_count; % have multiple paralogs compared with panmodel
GD_sort =  sum(GD_rxn ~= 0,1);% rxn number with GD based on speciessortclade

%% generalist/special list
for i = 1:length(strains_sortclade)
       cd(inputpath)
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    for j = 1:length(rxn)
      rxn_temp = split(rxn{j},';');
        [~,idx] = ismember(rxn_temp,model.rxns);
        genes = sum(model.rxnGeneMat(idx(idx~=0),:),1);
        [~,idx] = ismember(model.genes(find(genes)),OG_rxn_count(:,1));
        gene_counts = cell2mat(OG_rxn_count(idx(idx~=0),2));
        if ~isempty(idx==0) % deal with specific species genes.
            genes_temp = model.genes(find(genes));
            [~,id] = ismember(genes_temp(idx == 0),model.genes);
            temp = sum(model.rxnGeneMat(:,id),1);
            gene_counts(idx == 0) = temp;
        end
        general_rxn(j,i) = sum(gene_counts > 1); % gene number with muliple reactions 
    end
end

% Fig1a
subplot(1,6,1)
hold on
%f1 = scatter(x_total,clade_av,15,'filled','k');f1.MarkerFaceAlpha = 0.4;
h1 = boxplot(clade_av,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(h1,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('substrates','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
subplot(1,6,2)
%f2 = scatter(x_total,gain,15,'filled','k');f2.MarkerFaceAlpha = 0.4;
hold on
h2 = boxplot(gain,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[255,127,0]/255,'Labels',clades);
set(h2,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('gained','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
subplot(1,6,3)
hold on
%f3 = scatter(x_total,loss,15,'filled','k');f3.MarkerFaceAlpha = 0.4;
h3 = boxplot(loss,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[49,163,84]/255,'Labels',clades);
set(h3,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('lost ','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
subplot(1,6,4)
hold on
%f4 = scatter(x_total,HGT_species,15,'filled','k');f4.MarkerFaceAlpha = 0.4;
h4 = boxplot(HGT_species,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[215,48,39]/255,'Labels',clades);
set(h4,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Substrate HGT','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
subplot(1,6,5)
hold on
%f5 = scatter(x_total,GD_sort,15,'filled','k');f5.MarkerFaceAlpha = 0.4;
h5 = boxplot(GD_sort,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[118,42,131]/255,'Labels',clades);
set(h5,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Gene ratio','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
subplot(1,6,6)
hold on
%f6 = scatter(x_total,sum(general_rxn,1),15,'filled','k');f6.MarkerFaceAlpha = 0.4;
h6 = boxplot(sum(general_rxn,1),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[221,28,119]/255,'Labels',clades);
set(h6,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('General list','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
 %% Fig1d for HGT/GD/Generalist in metabolic gain
 for i = 1:length(strains_sortclade)
     gain_idx = split(gain_result(1,i),';');
     gain_idx = replace(gain_idx,SubBiologName(8:end),SubModelName(8:end));
     gain_idx = intersect(gain_idx,SubModelName);
     gain_idx(strcmp(gain_idx,'')) = [];
     HGT_idx = unique(trait(find(HGT_rxn(:,i)>=1)));
     GD_idx = unique(trait(find(GD_rxn(:,i)>=1)));
     GL_idx = unique(trait(general_rxn(:,i)>=1));
     contribution(i,1) = length(intersect(gain_idx,HGT_idx)); %  HGT 
     contribution(i,2) = length(intersect(gain_idx,GD_idx)); %  GD
     contribution(i,3) = length(intersect(gain_idx,GL_idx)); %  GL
     contribution(i,4) = length(setdiff(gain_idx,union(union(HGT_idx,GD_idx),GL_idx))); % no reason found
 end

 for i = 1:4
subplot(1,8,i*2-1)
h = boxplot(contribution(:,i),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
hold on
%f4 = scatter(x_total,contribution(:,i),'filled','k');f4.MarkerFaceAlpha = 0.4;
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('No. enzyme with paralog','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
end

% second figure
for i = 1:length(clades)
    idx = ismember(group,i);
    contribution_clade(i,1) = sum(contribution(idx,1),1);
    contribution_clade(i,2) = sum(contribution(idx,2),1);
    contribution_clade(i,3) = sum(contribution(idx,3),1);
    contribution_clade(i,4) = sum(contribution(idx,4),1);
end
h = bar(contribution_clade); %%%%%%%%%%%%%%%%%%%%%need to modify
set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Clade specific metabolic innovation','FontSize',12,'FontName','Helvetica');
legend({'HGT','Paralog','Generalist','Others'})
set(gca,'FontSize',10,'XTickLabelRotation',90)


% generate a figure of HGT/critical genes
 for i = 1:length(strains_sortclade)
     gain_idx = split(gain_result(1,i),';');
     gain_idx = replace(gain_idx,SubBiologName(8:end),SubModelName(8:end));
     gain_idx = intersect(gain_idx,SubModelName);
     gain_idx(strcmp(gain_idx,'')) = [];
     HGT_idx = unique(trait(find(HGT_rxn(:,i)>=1)));
     GD_idx = unique(trait(find(GD_rxn(:,i)>=1)));
     GL_idx = unique(trait(general_rxn(:,i)>=1));
     contribution(i,1) = length(intersect(gain_idx,HGT_idx)); %  HGT 
     contribution(i,2) = length(intersect(gain_idx,GD_idx)); %  GD
     contribution(i,3) = length(intersect(gain_idx,GL_idx)); %  GL
     contribution(i,4) = length(setdiff(gain_idx,union(union(HGT_idx,GD_idx),GL_idx))); % no reason found
 end

% Fig 1d
h1 = bar(sum(contribution_clade,1),'FaceColor',[255,127,0]/255,'FaceAlpha',0.3,'EdgeColor',[255,127,0]/255,'LineWidth',1);
set(gca,'XTick',1:1:4);
set(gca,'XTickLabel',{'HGT','Paralog','Generalist','Others'});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('No. subtrate gain','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);

% for all HGT where they locate
for i = 1:length(index(:,1)) % calculate which one is the first step
if endsWith(index{i},'_exp')
k = 0;
elseif startsWith(index{i},'r_')
    k = k+1;
    step = [step;k];
end
end
HGT_rxn_sum = sum(HGT_rxn,2);
sum(HGT_rxn_sum(step==1,:))

%% Figure for metabolic loss
% All HGT
fid2 = fopen('../data/accessory_new_HGT.tsv');
format = repmat('%s ',1,12);
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
all_HGT(:,i) = temp{i};
end
fclose(fid2);
all_HGT = all_HGT(find(strcmp(all_HGT(:,12),'Yes')),:);
temp = split(all_HGT(:,3),'@'); % split by the gene
all_HGT(:,3) = temp(:,1);

temp = split(all_HGT(:,7),'-');
all_HGT(:,12) = temp(:,4);
origin = {'Actinoba';'Amoebozoa';'Apicomple';'BacteroC';'Chlorofl';'Cyanobacte';'Euglenozo';'Firmicut';'Gemmatimo';'Metazoa';'Proteoba';'Streptoph';'Thaumarc';'other_Ba';'other_Stra'};
for j = 1:length(strains_sortclade)
    j
    ID = ismember(all_HGT(:,3),strains_sortclade(j,1));
    HGT_species(j) = length(unique(all_HGT(ID,2)));
    HGT_average_AI(j) = mean(cellfun(@str2num, all_HGT(ID,10)));
    for k = 1:length(origin)
        temp = find(startsWith(all_HGT(ID,12),origin{k}));
        origin_HGT(j,k) = length(temp);
    end
        
end
figure 
hold on
h2 = boxplot(HGT_species,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
f2 = scatter(x_total,HGT_species,20,'k','filled');f2.MarkerFaceAlpha = 0.4;
set(gca,'FontSize',10,'XTickLabelRotation',90)
ylabel('All metabolic HGTs in each clade','FontSize',12,'FontName','Helvetica');

%% Fig2b ALL ORIGIN OF hgt
temp = sum(origin_HGT,1);
[R C]=sort(temp,'descend');
lable = [origin(C(1:5));'Others'];
origin_HGT_final = [origin_HGT(:,C(1:5)),sum(origin_HGT(:,C(6:end)),2)];
for k = 1:6
subplot(1,6,k)
h = boxplot(origin_HGT_final(:,k),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(h,{'linew'},{1});
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
ylabel(strrep(lable{k},'_',' '),'FontSize',12,'FontName','Helvetica');
end

%% substrate loss
% define the metabolic loss as critical enzyme loss
[~,~,index] = xlsread('../data/substrateUsageGene1.xlsx','index');
[~,idx] = ismember(strains_sortclade,strainlist);
result_table = result_table(idx,:);
rxnexistence = result_table(:,startsWith(index(:,1),'r_'));
FBAresult = result_table(:,endsWith(index(:,1),'_model'));
sub_exp = result_table(:,endsWith(index(:,1),'_exp'));
index = index(startsWith(index(:,1),'r_'),:);

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

for i = 1:length(strains_sortclade) % define the loss type
    loss_temp = find(cell2mat(sub_exp(i,:))==0);
    for j = 1:length(loss_temp)
       rxn_idx = find(step(:,1)==loss_temp(j));
       criticalrxn_idx = cell2mat(index(step(:,1)==loss_temp(j),3)) > 0.7 & cell2mat(index(step(:,1)==loss_temp(j),4)) > 0.5;% based on the creteria of galactose
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
       criticalrxn_idx = cell2mat(index(step(:,1)==loss_temp(j),3)) > 0.7 & cell2mat(index(step(:,1)==loss_temp(j),4)) > 0.5;% based on the creteria of galactose
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
for i = 1:length(strains_sortclade) % define the loss type
    loss_temp = find(cell2mat(sub_exp(i,:))==0);
    for j = 1:length(loss_temp)
       rxn_idx = find(step(:,1)==loss_temp(j));
            temp = rxnexistence(i,rxn_idx); 
            model_temp = FBAresult(i,loss_temp(j));
            if cell2mat(model_temp) == 0 && all(cell2mat(temp))
                downsteam(i,loss_temp(j)) = 1;
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
     contribution(i,3) = length(intersect(loss_idx,dp_idx)); %  dp
     contribution(i,4) = length(intersect(loss_idx,general_idx)); %  general_list
     contribution(i,5) = length(setdiff(loss_idx,union(union(critical_idx,no_critical_idx),union(dp_idx,general_idx)))); % no reason found
 end

h1 = bar(sum(contribution,1),'FaceColor',[49,163,84]/255,'FaceAlpha',0.3,'EdgeColor',[49,163,84]/255,'LineWidth',1);
set(gca,'XTick',1:1:5);
set(gca,'XTickLabel',{'Critical rxn','Non-critical rxn','Downstream','Generalist','Others'});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('No. subtrate loss','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);
set(gca,'FontSize',10,'FontName','Helvetica');

clearvars general_idx dp_idx loss_idx no_critical_idx critical_idx temp rxnexistence result norxn_idx nocritical_noloss no_critical_noloss no_critical_loss loss_temp j k
clearvars reducedModel nocriticalrxn_idx m idx ID id i h h1 genes_temp genes sub_exp temp_no trait trait_1 trait_tmp trait_total format fid2 existence
%% Fig 2c critical enzyme
critical_idx = find(cell2mat(index(:,3)) > 0.7 & cell2mat(index(:,4)) > 0.5); % based on the creteria of galactose
% get trait
for i = 1:length(index(:,1))
if endsWith(index{i},'_exp')
trait_1 = replace(index{i,1},'_exp','');
trait{i} = trait_1;
else
trait{i} = trait_1;
end
end

h = bar([length(unique(trait(critical_idx))),32-length(unique(trait(critical_idx)))],'FaceColor','k','FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
set(gca,'FontSize',10,'FontName','Helvetica');
ylim([0 20])
xticks([1,2]);
xticklabels({'With critical reaction','Without critical reaction'})
ylabel('No. substrates','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(60);

%% Fig2f false positives
[~,~,index] = xlsread('../data/substrateUsageGene1.xlsx','index');
[~,idx] = ismember(strains_sortclade,strainlist);
result_table = result_table(idx,:); % sort the result table based on the strains_sortclade
rxnexistence = result_table(:,startsWith(index(:,1),'r_'));
FBAresult = result_table(:,endsWith(index(:,1),'_model'));
sub_exp = result_table(:,endsWith(index(:,1),'_exp'));

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
h = bar([sum(count_nomatch_generalist),sum(count_nomatch_nogeneralist)],'FaceColor','k','FaceAlpha',0.3,'EdgeColor','k','LineWidth',1);
ylabel('No. false positive','FontSize',10,'FontName','Helvetica','Color','k');
set(gca,'FontSize',10,'XTickLabelRotation',90)