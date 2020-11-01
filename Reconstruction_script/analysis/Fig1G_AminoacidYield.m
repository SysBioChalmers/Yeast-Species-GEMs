% This function is to generate  the figure of amino acid yield and ATP yield of different
% species.

%load panmodel
load('../ModelFiles/model_original.mat')
inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
Strain_information(:,i) = data{i};
end
fclose(fid2);
strains = Strain_information(:,1);

cd ../Reconstruction/otherchanges/
[~,~,mets_test,~,solresult] = getprecursorMatrixCobra(model_original,strains,inputpath);
yield_precursor = zeros(length(solresult),2);
save
% calculate the range of yield for each precursors
for i = 1:length(solresult)
    yield_precursor(i,1) = max(solresult(:,i));
    yield_precursor(i,2) = min(solresult(:,i));
end

% amino acids
mets = {'L-alanine [extracellular]';'L-arginine [extracellular]';'L-asparagine [extracellular]';'L-aspartate [extracellular]';...
    'L-cysteine [extracellular]';'L-glutamine [extracellular]';'L-glutamate [extracellular]';'L-glycine [extracellular]';...
    'L-histidine [extracellular]';'L-isoleucine [extracellular]';'L-leucine [extracellular]';'L-lysine [extracellular]';...
    'L-methionine [extracellular]';'L-phenylalanine [extracellular]';'L-proline [extracellular]';'L-serine [extracellular]';...
    'L-threonine [extracellular]';'L-tryptophan [extracellular]';'L-tyrosine [extracellular]';'L-valine [extracellular]'};
legend = {'L-ala';'L-arg';'L-asp';'L-asp';'L-cys';'L-glu';'L-glu';'L-gly';'L-his';'L-iso';'L-leu';'L-lys';'L-met';'L-phe';'L-pro';'L-ser';'L-thr';'L-try';'L-tyr';'L-val'};

[~,idx] = ismember(mets_test,model_original.mets);
metname = model_original.metNames(idx);

% get index of aa in biomass precursor
[~,idx] = ismember(mets,metname); 

% only aa
yield_aa = yield_precursor(idx,:);

% plot
spider_plot_R2019b(yield_aa',...
    'AxesInterval', 4,...
    'AxesDisplay', 'one',...
    'AxesPrecision', 1,...
    'AxesLimits', [zeros(1,20); repmat(3.2,1,20)],...
    'FillOption', 'on',...
    'FillTransparency', 0.1,...
    'LineWidth', 4,...
    'Marker', 'none',...
    'AxesFontSize', 14,...
    'AxesLabels', legend',...
    'LabelFontSize', 14);

%% find the ATP production by clade
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(data)
    Strain_information(:,i) = data{i};
end
fclose(fid2);
strainlist = Strain_information(:,1);
solresult = zeros(length(strainlist),1);
for i = 1:length(strainlist)
    cd(inputpath)
    load([strainlist{i},'.mat']);
    model = reducedModel;
    model= setParam(model,'ub',{'r_4046'},1000);
    model = setParam(model,'obj',{'r_4046'},1);
    sol = optimizeCbModel(model);
    solresult(i,1) = sol.f;
end

clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
result = [];
for i = 1:length(clades)
    idx = ismember(Strain_information(:,2),clades(i));
    result= [result;num2cell(solresult(idx,1)),repmat(clades(i),length(find(idx)),1)];
end

% box plot
h = boxplot(cell2mat(result(:,1)),result(:,2),'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255,'Labels',clades);
set(h,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
ylabel('ATP yielf on glucose','FontSize',10,'FontName','Helvetica','Color','k');
set(gcf,'position',[200 0 350 300]);
set(gca,'position',[0.11 0.31 0.77 0.65]);