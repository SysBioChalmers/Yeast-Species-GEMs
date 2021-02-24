%comparsion of different yeast models with publised ones

% This relies on anothe repo,please clone that one before use. https://github.com/SysBioChalmers/YeastsModels

inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';

% load mapping list published model/species
fid2 = fopen('../data/PublishedModel_list.tsv');
format = '%s %s %s %s %s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
data(:,i) = temp{i};
end
fclose(fid2);

current_path = pwd;
result = [];
result_model = [];
group = [];
for i = 1:length(data(:,1))
    result(1:3,i) = cellfun(@str2num,data(i,3:5));
    cd(inputpath)
    load([data{i,2},'.mat'])
    result(4:6,i) = [length(reducedModel.genes),length(reducedModel.rxns),length(reducedModel.mets)];
end
cd(current_path);

% plot
h = bar(result');
xlabel('Species','FontSize',10,'FontName','Helvetica','Color','k');
ylabel('Number','FontSize',10,'FontName','Helvetica','Color','k');
xticks(1:1:length(data(:,2)));
set(gca,'XTickLabel',strrep(data(:,2),'_',' '));
set(gca,'FontSize',10,'XTickLabelRotation',90)


