% This function is to plot whether the model can capture the evolution
% trait such as ura1 transfer and complexI loss

%% load data 
inputpath = '/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat/';
fid2 = fopen('../data/physiology/343_phenotype_clade.tsv');
format = '%s %s %s';
temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
for i = 1:length(temp)
Strain_information(:,i) = temp{i};
end
fclose(fid2);
strainlist = Strain_information(:,1);

%% ura1 gene transfer
current_path = pwd;
ura = zeros(343,3);
%ura('original','both','gene transfer')
for i = 1:343
    cd(inputpath)
    load([strainlist{i},'.mat']);
    model = reducedModel;
    [~,idx] = ismember('MNXR97425',model.rxnMetaNetXID); % ura1 transfer
    if idx~=0
        ura(i,3) = 1;
    end
    [~,idx] = ismember('MNXR97420',model.rxnMetaNetXID); %original ura9
    if idx ~=0
        ura(i,1) = 1;
    end
    if sum(ura(i,:))==2
        ura(i,2) = 1;
        ura(i,1) = 0;
        ura(i,3) = 0;
    end
end

clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};

for i = 1:length(clades)
idx = ismember(Strain_information(:,2),clades(i));
result(i,:)= sum(ura(idx,:),1);
end
color_set = [69,117,180
             252,141,89
             215,48,39]/255;
b = bar([result(:,1),result(:,2),result(:,3)],'stacked'); %ura1_original, ura1_both, ura1_m
    for k = 1:length(b)
        b(k).FaceColor = color_set(k,:);
        b(k).FaceAlpha = 0.8;
        b(k).EdgeColor = 'w';
        b(k).EdgeAlpha = 0;
    end

% mapping with crabtree effect
crabtree = {'Dekkera_bruxellensis';'Hanseniaspora_vinae';'Candida_glabrata';'Kluyveromyces_marxianus';'Lachancea_fermentati';...
    'Lachancea_kluyveri';'Lachancea_thermotolerans';'Lachancea_waltii';'Nakaseomyces_castellii';'Nakaseomyces_delphensis';'Naumovozyma_castellii';...
    'Saccharomyces_cerevisiae';'Saccharomyces_eubayanus';'Saccharomyces_kudriavzevii';'Saccharomyces_mikatae';'Saccharomyces_paradoxus';...
    'Saccharomyces_uvarum';'Tetrapisispora_phaffii';'Torulaspora_delbrueckii';'Vanderwaltozyma_polyspora';'Zygosaccharomyces_bailii';'Zygosaccharomyces_rouxii';...
    'yHMPu5000034876_Tetrapisispora_iriomotensis';'yHMPu5000034710_Kluyveromyces_dobzhanskii';'yHMPu5000026152_Torulaspora_franciscae';'Schizosaccharomyces_pombe'};
   [~,idx] = ismember(crabtree,Strain_information(:,1)); 
   [~,clade_idx] = ismember(Strain_information(idx,2),clades); % get all x
   for i = 1:length(idx)
       if ura(idx(i),1) ~=0
           y(i) = randi([0,floor(result(clade_idx(i),1)*100-1)],1,1)/100;
       elseif ura(idx(i),2) ~=0
           y(i) = randi([floor(result(clade_idx(i),1)*100 + 5),floor(sum(result(clade_idx(i),1:2))*100-5)],1,1)/100;
       elseif  ura(idx(i),3) ~=0
            y(i) = randi([floor(result(clade_idx(i),2)*100+1),100],1,1)/100;
       end
   end
hold on
scatter(clade_idx,y,50,'x','k')
set(gca,'FontSize',10,'FontName','Helvetica');
% set(gca,'ycolor','k');
ylabel('URA1 source analysis','FontSize',10,'FontName','Helvetica','Color','k');
set(gca,'XTickLabel',clades);
set(gca,'FontSize',10,'XTickLabelRotation',90)

% map anaerobic condition with 
anaerobic = {'Sugiyamaella_lignohabitans';'Dekkera_bruxellensis';'yHMPu5000034625_Pichia_kudriavzevii';'yHMPu5000026142_Citeromyces_matritensis';'Candida_albicans';'Candida_parapsilosis';'Candida_tropicalis';...
             'Clavispora_lusitaniae';'Spathaspora_passalidarum';'Wickerhamia_fluorescens';'Wickerhamomyces_anomalus';'yHMPu5000035686_Cyberlindnera_saturnus';'Hanseniaspora_uvarum';'Hanseniaspora_valbyensis';...
             'Hanseniaspora_vinae';'yHMPu5000034957_Hanseniaspora_osmophila';'Ashbya_aceri';'Candida_glabrata';'Eremothecium_coryli';'Kluyveromyces_lactis';'Kluyveromyces_marxianus';'Lachancea_fermentati';...
             'Lachancea_kluyveri';'Lachancea_thermotolerans';'Lachancea_waltii';'Nakaseomyces_bacillisporus';'Nakaseomyces_castellii';'Nakaseomyces_delphensis';'Naumovozyma_castellii';'Naumovozyma_dairenensis';...
             'Saccharomyces_cerevisiae';'Saccharomyces_eubayanus';'Saccharomyces_paradoxus';'Saccharomyces_uvarum';'Tetrapisispora_blattae';'Tetrapisispora_phaffii';'Torulaspora_delbrueckii';'Vanderwaltozyma_polyspora';...
             'Zygosaccharomyces_bailii';'yHAB154_Kazachstania_transvaalensis';'yHMPu5000034881_Torulaspora_pretoriensis';'yHMPu5000034876_Tetrapisispora_iriomotensis';'yHMPu5000034862_Zygotorulaspora_florentina';...
             'yHMPu5000026152_Torulaspora_franciscae';'Schizosaccharomyces_pombe'};

         [~,idx] = ismember(anaerobic,Strain_information(:,1));
         [~,clade_idx] = ismember(Strain_information(idx,2),clades); % get all x
         for i = 1:length(idx)
             if ura(idx(i),1) ~=0
                 y(i) = randi([0,floor(result(clade_idx(i),1)*100)],1,1)/100;
             elseif ura(idx(i),2) ~=0
                 y(i) = randi([floor(result(clade_idx(i),1)*100),floor(sum(result(clade_idx(i),1:2))*100)],1,1)/100;
             elseif  ura(idx(i),3) ~=0
                 y(i) = randi([floor(result(clade_idx(i),2)*100),100],1,1)/100;
             end
         end

hold on
scatter(clade_idx,y,50,'x','k')
legend({'URA1[mito]','URA1[both]','URA1[c]','Crabtree postive species','Anaerobic species'},'FontSize',6,'FontName','Helvetica','location','se');

%% complex I existence
fid2 = fopen('../data/physiology/complexI_existence.tsv');
format = '%s %s';
data = textscan(fid2,format,'Delimiter','\t','HeaderLines',1);
strain_withoutcomplexI = data{1};
count = cellfun(@str2num,data{2});

[~,idx] = ismember(lower(strain_withoutcomplexI),lower(Strain_information(:,1)));
strain_withoutcomplexI(idx,:) = strain_withoutcomplexI;
count(idx,:) = count;
clades = {'Ascomycota';'Lipomycetaceae';'Trigonopsidaceae';'Dipodascaceae/Trichomonascaceae';'Alloascoideaceae';'Sporopachydermia';'Pichiaceae';'CUG-Ala';'CUG-Ser1';'CUG-Ser2';'Phaffomycetaceae';'Saccharomycodaceae';'Saccharomycetaceae'};
clades(1) = [];% take away the outgroup
result = [];
group = [];
for i = 1:length(clades)
idx = ismember(Strain_information(:,2),clades(i));
result = [result;count(idx)];
group = [group;repmat(i,length(count(idx)),1)];
end

% figure
h = boxplot(result,group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
set(h,{'linew'},{1});
ylabel('NADH:ubiquinone oxidoreductase subunit number','FontSize',10,'FontName','Helvetica','Color','k');
camroll(-90)