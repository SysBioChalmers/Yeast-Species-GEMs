currentpath = pwd;
%Figure S1a
fid2 = fopen('../data/genome_size_gene.tsv');
format = '%s %s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
Strain_information(:,i) = data{i};
end
fclose(fid2);

clades = unique(Strain_information(:,2));
numrxn = zeros(length(clades),1);
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
result = [];

[~,idx] = ismember(Strain_information(:,4),clades);
Strain_information(:,5) = num2cell(idx);

h = boxplot(cellfun(@str2num, Strain_information(:,3)),cell2mat(Strain_information(:,5)),'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');

set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gcf,'position',[200 0 350 300]);
set(gca,'position',[0.11 0.31 0.77 0.65]);

ylabel('Genome size','FontSize',14,'FontName','Helvetica','Color','k');

%% genrate the plot for rxn number in all species, not included in figures of paper
cd ../Reconstruction/otherchanges/
inputpath = '../modelRelated/ssGEMs';
load('../modelRelated/panModel.mat');
[~,rxnMatrix,~] = getprecursorMatrixCobra(model_original,Strain_information(:,1),inputpath,{''},0);
x = sum(rxnMatrix,1);
h = histogram(x,20,'FaceColor',[56,108,176]/255,'FaceAlpha',0.3);
set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
set(gcf,'position',[200 0 350 300]);
xlabel('Occurence number','FontSize',14,'FontName','Helvetica','Color','k');
ylabel('Reaction number','FontSize',14,'FontName','Helvetica','Color','k');
set(gca,'position',[0.21 0.31 0.77 0.65]);

%% figure for unique rxn without no comp
% unique metabolite to mets without comps
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

% get the count for each gene
result = [model_original.genes,num2cell(sum(newrxngenemat,1))'];
result_sorted = sort(cell2mat(result(:,2)));

xlim([0 50])
h = cdfplot(result_sorted);
set( h, 'linewidth',3);
set(gca,'FontSize',20,'FontName','Helvetica');
ylabel('Percentage','FontSize',24,'FontName','Helvetica','Color','k');
xlabel('OG linked reaction number','FontSize',24,'FontName','Helvetica','Color','k');
cd(currentpath)