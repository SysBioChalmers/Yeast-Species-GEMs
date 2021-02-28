% Figure 2b: this function is to generate plot for pan/core/accessory reactions Figure 2b and
% also the dataset for pathway figure using iPath3.0
currentpath= pwd;
load('../Reconstruction/modelRelated/StrainData.mat')
inputpath = '../Reconstruction/modelRelated/ssGEMs';
strains = StrianData.strains;
cd(inputpath)
for j = 1:10 % sample 10 times
    j
    samplesize = [1:1:343];
    for i = 1:343
        sample=samplesize(randperm(length(samplesize),1)); % random sample one strain
        samplesize = setdiff(samplesize,sample); % exclude the sampled one
        m = strains{sample};
        load([m,'.mat'])
        model = reducedModel;
        if i == 1
            core_rxn = model.rxns;
            acce_rxn = model.rxns;
            pan_rxn = model.rxns;
            core_final(i,j) = length(core_rxn);
            acce_final(i,j) = length(acce_rxn);
            pan_final(i,j) = length(pan_rxn);
        else
            core_rxn = intersect(core_rxn,model.rxns);
            pan_rxn = union(model.rxns,pan_rxn);
            acce_rxn = setdiff(pan_rxn,core_rxn);
            core_final(i,j) = length(core_rxn);
            acce_final(i,j) = length(acce_rxn);
            pan_final(i,j) = length(pan_rxn);
        end
    end
end

x1 = mean(acce_final,2); x1(1) = 0;
x2 = mean(core_final,2);
x3 = mean(pan_final,2);
figure
hold on
plot([1:1:343],x1,'k-','LineWidth',2,'color',[215,48,39]/255)
plot([1:1:343],x2,'k-','LineWidth',2,'color','k')
plot([1:1:343],x3,'k-','LineWidth',2,'color',[69,117,180]/255)
set(gca,'FontSize',6,'FontName','Helvetica');
set(gca,'ycolor','k');
set(gca,'xcolor','k');
ylabel('Reaction number','FontSize',8,'FontName','Helvetica','Color','k');
xlabel('Sampled species','FontSize',8,'FontName','Helvetica','Color','k');
legend({'Accessory','Core','Pan'},'FontSize',6,'FontName','Helvetica','location','se');
xlim([0,350])
set(gcf,'position',[10 50 200 100]);
set(gca,'position',[0.15 0.2 0.65 0.85]);

% generate data for ipath3.0 for showing those rxns in pathway map
[~,ID] = ismember(acce_rxn,model_original.rxns);
acce_rxn_kegg = model_original.rxnKEGGID(ID(ID~=0));
[~,ID] = ismember(core_rxn,model_original.rxns);
core_rxn_kegg = model_original.rxnKEGGID(ID(ID~=0));
cd(currentpath)