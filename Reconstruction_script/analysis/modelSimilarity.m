
% this function is to get the model similarity analysis for all models we
% generated



% load group information

[a, b, Strain_information]=xlsread('../../ComplementaryData/genome_summary_332_yeasts.xlsx','clades');
Strain_information = Strain_information(2:end,:);
clades = unique(Strain_information(:,2));
for i = 1:length(strains)
    for j = 1:i
        tempa = model.rxns(logical(rxnMatrix(i,:)));
        tempb = model.rxns(logical(rxnMatrix(j,:)));
        jaccard_sim = length(intersect(tempa,tempb))/length(union(tempa,tempb));
        d(j,i) = 1 - jaccard_sim;
        d(i,j) = 1 - jaccard_sim;
    end
end

% we take sce as the first one and find all models toghther

sce = find(contains(strains,'Saccharomyces_cerevisiae'));
dis = d(sce,:);
[Y,I] = sort(dis);
[~,idx] = ismember(Strain_information(I,2),clades)

color_palette = [166,206,227
31,120,180
178,223,138
51,160,44
251,154,153
227,26,28
253,191,111
255,127,0
202,178,214
106,61,154
255,255,53; 150,33,62;1,5,6]/255;


for i = 1:length(I)
m = strains{I(i)};
load([m,'.mat'])
reducedModel.id = m;
models{i} = reducedModel;
end
compStruct = compareMultipleModels(models,true,true,Strain_information(I,2))

[~,idx] = ismember(Strain_information(I,2),clades);
b = bar(repmat(1,343,1),1);
b.FaceColor = 'flat';
for i = 1:13
x = find(idx == i);
b.CData(x,:) = repmat(color_palette(i,:),length(x),1);
end