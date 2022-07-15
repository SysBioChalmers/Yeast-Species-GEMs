function [transall,transselected] = sortTransGapseq(strains,inputpath,mets)
% This function is to sort transport annotation from gapseq for 332 yeast
% species; we would like to find whether transporters existing for those
% inputted mets
% inputpath = '../Reconstruction/gapseqresult'; % downloaded from figshare 
% load('../Reconstruction/modelRelated/StrainData.mat');
% strains = StrianData.strains;
% [~,~,index] = xlsread('../Reconstruction/modelRelated/substrateUsageGene.xlsx','index'); % this file is missing in current folder.
% mets = strrep(index(endsWith(index(:,1),'_exp')),'_exp','')
file = dir(inputpath);
file = {file.name};
file = setdiff(file,'.');
file = setdiff(file,'..');
transall = [];
for i = 1:length(file)
    disp([file{i},': ', num2str(i),'/',num2str(length(file))])
    fileName    = file{i};
    fID         = fopen([inputpath,'/',fileName]);
    tmp    = textscan(fID,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    fclose(fID);
    for j = 1:length(tmp)
        trans(:,j) = tmp{j};
    end
    trans(:,11) = strrep(trans(:,11),' ','&');
    trans(:,11) = extractAfter(trans(:,11),'gene=');
    strain = extractBefore(fileName,'.max-Transporter');
    strain = regexprep(strain,'_16(\d*)','');
    
    % transfer id into the strain@seq ID
    [~,idx_strain] = ismember(lower(strain),lower(strains));
    fileName    = ['../Multi_scale_evolution/pan_genome/result/id_mapping/',strains{idx_strain},'.tsv'];
    fID         = fopen(fileName);
    protData    = textscan(fID,'%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    fclose(fID);
    
    geneID_core = protData{2}; % geneID in 343 yeast species with @seq
    panID_final = protData{5}; % panID
    OG_ID = protData{4}; 
    draftgeneID      = protData{7}; % mRNAID geneID CDS
    draftgeneID      = strrep(draftgeneID,' ','&'); % mRNAID&geneID&CDS
    draftgeneID = extractAfter(draftgeneID,'gene=');
    panID_final = strrep(panID_final,'Saccharomyces_cerevisiae@','');
    
    [~,idx] = ismember(trans(:,11),draftgeneID);
    trans(:,14) = strains(idx_strain);
    trans(:,15) = geneID_core(idx);
    trans(:,16) = panID_final(idx);
    trans(:,17) = OG_ID(idx);
    transall = [transall;trans];
    clear trans;
end

% 
mets = mets(~cellfun(@isempty,mets(:,2)),:);
[~,idx] = ismember(transall(:,3),mets(:,2));
transselected = transall(idx~=0,:);

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
strains_sortclade = [];

for i = 1:length(clades)
idx = ismember(Strain_information(:,2),clades(i));
[~,ID] = ismember(Strain_information(idx,1),strains);
group = [group;i*ones(length(ID(ID~=0)),1)];
strains_sortclade = [strains_sortclade;strains(ID(ID~=0))];
end

a = tabulate(transselected(:,14));
b = tabulate(transall(:,14));
[~,ID] = ismember(strains_sortclade,a(:,1));
[~,ID2] = ismember(strains_sortclade,b(:,1));
clade_av = cell2mat(a(ID,2)); % trans according to clade
clade_av2 = cell2mat(b(ID2,2)); % trans according to clade
subplot(1,2,1)
h3 = boxplot(clade_av,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[49,163,84]/255,'Labels',clades);
set(h3,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Selected Transport rxns','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',0)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
hold on
camroll(-90)
subplot(1,2,2)
h4 = boxplot(clade_av2,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[215,48,39]/255,'Labels',clades);
set(h4,{'linew'},{1});
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('Transport rxns','FontSize',12,'FontName','Helvetica');
set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
camroll(-90)
end
