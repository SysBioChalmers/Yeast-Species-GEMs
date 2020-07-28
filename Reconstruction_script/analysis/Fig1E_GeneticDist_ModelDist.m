% This function is to generate the figure 1e

% load group information
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
    Strain_information(:,i) = data{i};
end
fclose(fid2);
clades = unique(Strain_information(:,2));
clear data
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
% we take sce as the first one and find all models toghther%
strain_gene = genetic(1,2:end);
genetic = genetic(2:end,2:end);

sce = find(contains(strain_gene,'yHMPu5000034963_Hanseniaspora_clermontiae'));
%sce = find(contains(strain_gene,'Schizosaccharomyces_pombe'));
dis_gene = cellfun(@str2num, genetic(sce,1:end));

% find the similarity of the phenotype data
fid2 = fopen('../data/physiology/Biolog_substrate.tsv');
format = repmat('%s ',1,333);
format = strtrim(format);
substrate = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(substrate)
    data(:,i) = substrate{i};
end
strainlist = data(1,5:end);

sce = find(contains(strainlist,'yHMPu5000034963_Hanseniaspora_clermontiae'));
data(1,:) = [];
SubModelName = data(:,2);
data(:,1:4) = [];
fclose(fid2);

for i = 1:length(strainlist)
    for j = 1:length(strainlist)
    tempa = data(1:59,i);
    tempb = data(1:59,j);
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
    temp = 1-pdist([tempa,tempb]','jaccard');
    dis_p(i,j) = temp(1,1);
    dis_p(j,i) = temp(1,1);
    end
end
sce = find(contains(strainlist,'yHMPu5000034963_Hanseniaspora_clermontiae'));
%sce = find(contains(strainlist,'Schizosaccharomyces_pombe'));

% figure out the model predicted phenotype distance % FBAresult is based on
% strainlist
load('../Reconstruction/FBAresult6.mat')
FBAresult(find(abs(FBAresult)>0)) = 1;
for i = 1:length(strainlist)
    for j = 1:length(strainlist)
    tempa = FBAresult(1:59,j);
    tempb = FBAresult(1:59,i);
    noexista = find(isnan(tempa));
    noexistb = find(isnan(tempb));
    noexist = union(noexista,noexistb);
    tempa(noexist) = [];
    tempb(noexist) = [];
    temp = 1-pdist([tempa,tempb]','jaccard');
    dis_p_model(i,j) = temp(1,1);
    dis_p_model(j,i) = temp(1,1);
    end
end

% compare the distance with sce and other strains
[~,ID] = ismember(strainlist,strain_gene);
dis_gene_new = dis_gene(ID); % reorder the dis_rxn based on the phenotype

%plot https://blog.csdn.net/weixin_42943114/article/details/90074259
%??????
Nx=500;
Ny=500;
X = dis_gene_new;
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
