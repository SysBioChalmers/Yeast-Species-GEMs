% This function is to generate the figure 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2B clade specific reactions & lpot the wind rose figure according
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
sub = replace(sub,'Gluconeogenesis;','')
sub = replace(sub,'Metabolic pathways;','')


    
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
[a,b] = sort(sum(matrix,1),'descend')
subID = subID(b)
matrix = matrix(:,b)


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure rxn comparsion for each clade
clades = unique(Strain_information(:,2));
group = [];
clade_rxn = [];
rxnnum = sum(rxnMatrix,2);
strainlist = StrianData.strains;
for i = 1:length(clades)
idx = ismember(Strain_information(:,2),clades(i));
[~,ID] = ismember(Strain_information(idx,1),strainlist);
clade_rxn = [clade_rxn;rxnnum(ID(ID~=0))];
group = [group;(i-1)*ones(length(rxnnum(ID(ID~=0))),1)];
end
h = boxplot(clade_rxn,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');

set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gcf,'position',[200 0 350 300]);
set(gca,'position',[0.11 0.31 0.77 0.65]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2C Substrate Usage test
% load group information
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
Strain_information(:,i) = data{i};
end
fclose(fid2);
clades = unique(Strain_information(:,2));

%load models
current_path = pwd;
% cd('../ModelFiles/mat')
% for i = 1:length(strains)
%     m = strains{i};
%     load([m,'.mat'])
%     reducedModel.id = m;
%     models{i} = reducedModel;
% end
% cd(current_path)
% 
% % Create binary matrix of reactions
% [id,binary_matrix] = compareModelField(models,'rxns');
% compStruct.reactions.IDs = id;
% compStruct.reactions.matrix = binary_matrix;
% 
% % calculate jaccard similarity
% compStruct.structComp = 1-squareform(pdist(binary_matrix','jaccard'));
% 
% % we take sce as the first one and find all models toghther
% sce = find(contains(strains,'Saccharomyces_cerevisiae'));
% sce = find(contains(strains,'Schizosaccharomyces_pombe'));
% dis_rxn = compStrucxt.structComp(sce,:);

% load genetic distance
fid2 = fopen('../data/taxa_pairwise_dist_343.txt');
format = repmat('%s ',1,344);
format = strtrim(format);
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(temp)
genetic(:,i) = temp{i};
end

% we take sce as the first one and find all models toghther
%sce = find(contains(genetic(1,2:end),'Saccharomyces_cerevisiae'));
sce = find(contains(genetic(1,2:end),'Schizosaccharomyces_pombe'));
dis_gene = cellfun(@str2num, genetic(sce,2:end));



% find the similarity of the phenotype data
fid2 = fopen('../data/physiology/Biolog_substrate.tsv');
format = repmat('%s ',1,333);
format = strtrim(format);
substrate = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(substrate)
    data(:,i) = substrate{i};
end
strainlist = data(1,5:end);
data(1,:) = [];
SubModelName = data(:,2);
data(:,1:4) = [];
fclose(fid2);


for i = 1:length(strainlist)
    tempa = data(1:59,330);
    tempb = data(1:59,i);
    tempa = strrep(tempa,'n','NAN');
    tempb = strrep(tempb,'n','NAN');
    tempa = strrep(tempa,'v','1');
    tempb = strrep(tempb,'v','1');
    noexista = find(strcmp(tempa,'NAN'));
    noexistb = find(strcmp(tempb,'NAN'));
    noexist = union(noexista,noexistb);
    tempa(noexist) = [];
    tempb(noexist) = [];
    tempa = cellfun(@str2num, tempa);
    tempb = cellfun(@str2num, tempb);
    temp = 1-squareform(pdist([tempa,tempb]','jaccard'));
    dis_p(i,1) = temp(1,2);
end

%sce = find(contains(strainlist,'Saccharomyces_cerevisiae'));
sce = find(contains(strainlist,'Schizosaccharomyces_pombe'));

% figure out the model predicted phenotype distance
for i = 1:length(strainlist)
    tempa = FBAresult(1:59,sce);
    tempb = FBAresult(1:59,i);
    noexista = find(strcmp(tempa,'n'));
    noexistb = find(strcmp(tempb,'n'));
    noexist = union(noexista,noexistb);
    tempa(noexist) = [];
    tempb(noexist) = [];
    tempa = cellfun(@str2num, tempa);
    tempb = cellfun(@str2num, tempb);
    temp = 1-squareform(pdist([tempa,tempb]','jaccard'));
    dis_p_model(i,1) = temp(1,2);
end

% compare the distance with sce and other strains
[~,ID] = ismember(strainlist,genetic(2:end,1));
dis_gene = dis_gene(ID); % reorder the dis_rxn based on the phenotype



%plot https://blog.csdn.net/weixin_42943114/article/details/90074259
%??????
Nx=500;
Ny=500;
X = dis_gene;
Y = dis_p';
hold on 
Xedge=linspace(min(X),max(X),Nx);
Yedge=linspace(min(Y),max(Y),Ny);

%N?xy??????
[N,~,~,binX,binY] = histcounts2(X,Y,[-inf,Xedge(2:end-1),inf],[-inf,Yedge(2:end-1),inf]);

XedgeM=movsum(Xedge,2)/2;
YedgeM=movsum(Yedge,2)/2;

[Xedgemesh,Yedgemesh]=meshgrid(XedgeM(2:end),YedgeM(2:end));

%??pcolor?
figure(1)
pcolor(Xedgemesh,Yedgemesh,N');shading interp

%????
%h=ones(round(Nx/20));
%h=fspecial('disk',round(Nx/40));
h = fspecial('gaussian',round(Nx/20),6);%????????
N2=imfilter(N,h);
figure(2)
pcolor(Xedgemesh,Yedgemesh,N2');shading interp

ind = sub2ind(size(N2),binX,binY);
col = N2(ind);

figure(3)
scatter(X,Y,50,col,'filled','MarkerFaceAlpha',0.7);
ylabel('Phenotypic similarity','FontSize',12,'FontName','Helvetica','Color','k');
xlabel('Genetic distance','FontSize',12,'FontName','Helvetica','Color','k');
set(gcf,'position',[0 0 170 190]);
set(gca,'position',[0.19 0.1 0.78 0.85]);
P = polyfit(X,log(Y),1);
yi = exp(polyval(P,linspace(min(X),max(X),100)));

hold on 
plot(linspace(min(X),max(X),100),yi)
scatter(X,dis_p_model',50,'s','filled','MarkerFaceAlpha',0.7)
P = polyfit(X,log(dis_p_model'),1);
yi = exp(polyval(P,linspace(min(X),max(X),100)));
plot(linspace(min(X),max(X),100),yi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2E
[a, b, Strain_information]=xlsread('../../ComplementaryData/genome_summary_332_yeasts.xlsx','clades');
Strain_information = Strain_information(2:end,:);
clades = unique(Strain_information(:,2));
strains = Strain_information(:,1);
compStruct = compareMultipleModels(models,true,true,Strain_information(I,2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panrxn core rxn accesory rxn

for i = 1:343
m = StrianData.strains{i};
load([m,'.mat'])
model = reducedModel;
if i == 1
core_rxn = model.rxns;
acce_rxn = model.rxns;
pan_rxn = model.rxns;
a(i,1) = length(core_rxn);
a(i,2) = length(acce_rxn);
a(i,3) = length(pan_rxn);
else
core_rxn = intersect(core_rxn,model.rxns);
pan_rxn = union(model.rxns,pan_rxn);
acce_rxn = setdiff(pan_rxn,core_rxn);
a(i,1) = length(core_rxn);
a(i,2) = length(acce_rxn);
a(i,3) = length(pan_rxn);
end
end

hold on
plot([1:1:343],a(:,3),'k-','LineWidth',2,'color',[215,48,39]/255)
plot([1:1:343],a(:,1),'k-','LineWidth',2,'color','k')
plot([1:1:343],a(:,2),'k-','LineWidth',2,'color',[69,117,180]/255)
set(gca,'FontSize',12,'FontName','Helvetica');
ylabel('Reaction number','FontSize',20,'FontName','Helvetica','Color','k');
xlabel('Sampled speceis model','FontSize',20,'FontName','Helvetica','Color','k');
xlabel('Sampled speceis models','FontSize',20,'FontName','Helvetica','Color','k');
set(gca,'FontSize',18,'FontName','Helvetica');
ylabel('Reaction number','FontSize',20,'FontName','Helvetica','Color','k');
xlabel('Sampled speceis models','FontSize',20,'FontName','Helvetica','Color','k');

%%%%%%%%
% Addtional function
function [id,compMat] = compareModelField(models,field)
    % Generates a list of unique field entries and a matrix quantifying the
    % number of appearances of each field entry in each model
    
    % get unique list of field entries
    hasfield = cellfun(@(m) isfield(m,field),models);
    id = cellfun(@(m) m.(field),models(hasfield),'UniformOutput',false);
    id = unique(vertcat(id{:}));
    
    % assemble matrix comparing frequency of each entry in each model
    compMat = zeros(numel(id),numel(models));
    for i = 1:numel(models)
        [~,entryIndex] = ismember(models{i}.(field),id);  % get index of each field entry in the unique id list
        compMat(:,i) = histcounts(entryIndex, 0.5:1:(numel(id)+0.5));  % determine the frequency at which each index appears
    end
end