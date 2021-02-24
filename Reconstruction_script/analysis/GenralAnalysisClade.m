% Figure 2a: This function is to generate figure 1 for analysis of models on species number, reaction number, accessory reaction
% number, substrate usage number and biomass yield based on clades
currentpath = pwd;
% Figure 1
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
    Strain_information(:,i) = temp{i};
end
fclose(fid2);

clades = unique(Strain_information(:,2));
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
group = [];
clade_av = [];

for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    clade_species(i) = length(find(idx)); % species number in each clade
end

%% for species based on clade
subplot(1,6,1)
h = bar(clade_species,'LineWidth',1);
h.FaceColor = [56,108,176]/255;
h.EdgeColor = [56,108,176]/255;
h.FaceAlpha = 0.3;
h.BarWidth = 0.5;
ylabel('No. species','FontSize',10,'FontName','Helvetica','Color','k');
set(gca,'xtick',[])
camroll(-90)

%%  for reaction number based on clade
cd ../Reconstruction/otherchanges
strains_sortclade = [];
for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    group = [group;i*ones(length(idx(idx~=0)),1)];
    strains_sortclade = [strains_sortclade;Strain_information(idx~=0,1)];
end

load('../modelRelated/panModel.mat');
inputpath = '../modelRelated/ssGEMs';
[~,rxnMatrix,~,~,~] = getprecursorMatrixCobra(model_original,strains_sortclade,inputpath,[],0);
rxn_count = sum(rxnMatrix,2);

% figure for rxns based on clade
subplot(1,6,2)
h = boxplot(rxn_count,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
set(h,{'linew'},{1});
set(gca,'xtick',[])
ylabel('No. reactions','FontSize',10,'FontName','Helvetica','Color','k');
camroll(-90)

%%  for access reactions based on clade
rxn_accessory = rxnMatrix;
rxn_accessory(:,sum(rxn_accessory,1)> 343*0.95)=[];
acce_rxn_count = sum(rxn_accessory,2);
subplot(1,6,3)
h = boxplot(acce_rxn_count,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[255,127,0]/255);
set(h,{'linew'},{1});
set(gca,'xtick',[])
ylabel('No. accessory reactions','FontSize',10,'FontName','Helvetica','Color','k');
camroll(-90)

%%  for genes based on clade
for i = 1:length(strains_sortclade)
    cd(inputpath)
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    genes(i) = length(model.genes);
end
cd(current_path)
subplot(1,6,4)
h = boxplot(genes,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[49,163,84]/255,'Labels',clades);
set(h,{'linew'},{1});
set(gca,'xtick',[])
ylim([500,2000])
ylabel('No. genes','FontSize',10,'FontName','Helvetica','Color','k');
camroll(-90)

%% figure for number of substrates that can be ultilized 
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
clade_av = [];
for i = 1:length(strains_sortclade)
    [~,ID] = ismember(strains_sortclade(i,1),strainlist);
    data_slected = data(1:58,ID(ID~=0));
    if isempty(data_slected)
        count = zeros(1,1);
    else
        count = sum(cell2mat(data_slected),1,'omitnan');
    end
    clade_av = [clade_av;count'];
end
subplot(1,6,5)
h = boxplot(clade_av,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('substrates','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)

% figure 1f for substrate prediction based on clade
%figure is plotted in substrateanalysis

%% figure 1f for biomass yield
for i = 1:length(strains_sortclade)
    cd(inputpath)
    load([strains_sortclade{i},'.mat']);
    model = reducedModel;
    sol = solveLP(model);
    solresult(i) = sol.f;
end
subplot(1,6,6)
h = boxplot(abs(solresult),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[118,42,131]/255,'Labels',clades);
set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Biomass yield','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
cd(currentpath)