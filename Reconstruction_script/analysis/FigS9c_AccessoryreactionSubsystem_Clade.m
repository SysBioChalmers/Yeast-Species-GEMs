%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1e clade specific reactions & plot the wind rose figure according
% to the evolution clade

% exception code for new rxns subsystems
fileName = '/Users/feiranl/Box Sync/hongzhong/new_rxn_mnx_keggpathway.txt';
fID      = fopen(fileName);
mapping  = textscan(fID,'%s%s','Delimiter','\t','HeaderLines',1);
keggid   = mapping{1};
pathway  = mapping{2};

for i = 1:length(model_original.rxns)
    if strcmp(model_original.subSystems{i},'')
        [~,ID] = ismember(model_original.rxnKEGGID(i),keggid);
        if ID ~= 0
            model_original.subSystems{i} = pathway(ID);
        end
    end
end

%change subsystem into one cell
for i = 1:length(model_original.rxns)
    aaa = model_original.subSystems{i};
    bbb = join(aaa,';');
    sub(i,1) = bbb;
end
sub = replace(sub,'Gluconeogenesis;','');
sub = replace(sub,'Metabolic pathways;','');

fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
    Strain_information(:,i) = data{i};
end
fclose(fid2);
clades = unique(Strain_information(:,2));
numrxn = zeros(length(clades),1);
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
result = [];
for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    rxns_clade = rxnMatrix(idx,:);
    rxns_outclade_common = rxnMatrix(~idx,:);
    rxns_clade = sum(rxns_clade,1);
    rxns_clade = find(rxns_clade==length(find(idx))); % find union of those rxns if only common rxns use  == length(find(idx))
    corerxns = find(sum(rxnMatrix,1) == 343);
    %rxns_outclade_common = sum(rxns_outclade_common,1);
    %rxns_outclade_common = find(rxns_outclade_common>0); % find union of those rxns
    %unique_rxns = setdiff(rxns_clade,rxns_outclade_common);
    unique_rxns = setdiff(rxns_clade,corerxns);
    numrxn(i) = length(unique_rxns);
    subs_unique = join(sub(unique_rxns,1),';');
    subs_unique = split(subs_unique,';');
    order = tabulate(subs_unique);
    order_idx = sort(cell2mat(order(:,2)),'descend');
    idx = find(ismember(cell2mat(order(:,2)),order_idx(order_idx>2)));
    if ~isempty(idx)
        result= [result;order(idx,1:2),repmat(clades(i),length(idx),1)];
    end
end

% sort the result
subID = unique(result(:,1));
matrix = zeros(length(clades),length(subID));
for i = 1:length(subID)
    for j = 1:length(clades)
        mapping = find(strcmp(result(:,1),subID(i)) & strcmp(result(:,3),clades(j)));
        if ~isempty(mapping)
            matrix(j,i) = cell2mat(result(mapping,2));
        end
    end
end

% sort the data
[~,b] = sort(sum(matrix,1),'descend');
subID = subID(b);
matrix = matrix(:,b);

% plot the result
% will replace the bar plot with windrose plot
hold on
color_bg = [128,128,128]/255;
color = [165,0,38
    215,48,39
    244,109,67
    253,174,97
    254,224,144
    224,243,248
    171,217,233
    116,173,209
    69,117,180
    49,54,149]/255;
%h1 = plot([1:1:length(clades)],numrxn,'ko','LineWidth',1,'MarkerSize',5);
h1 = bar(numrxn,'FaceColor',color_bg,'FaceAlpha',0.3,'EdgeColor',color_bg,'LineWidth',0.5);
h2 = bar(matrix(:,1:10),'stacked','EdgeColor','k','LineWidth',0.5);
for i = 1:10
    h2(i).FaceAlpha = 0.6;
    h2(i).FaceColor = color(i,:);
end
set(gca,'XTick',1:1:length(clades));
set(gca,'XTickLabel',clades);
set(gca,'XTickLabelRotation',90)
expression = 'sce(\w+)\s';
replace = '';
subID = regexprep(subID,expression,replace);
subID = strtrim(subID);
legend([{'Others'},subID(1:10)'])
set(gca,'FontSize',8,'FontName','Helvetica');
set(gca,'ycolor','k');
set(gca,'xcolor','k');
ylabel('Clade shared accessory reaction','FontSize',10,'FontName','Helvetica','Color','k');
set(gcf,'position',[10 0 330 330]);
set(gca,'position',[0.15 0.33 0.65 0.55]);
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%