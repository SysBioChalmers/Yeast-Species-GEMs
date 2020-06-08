% this function is to generate plot for pan/core/accessory reactions and
% also the dataset for pathway figure using iPath3.0

load('../Reconstruction/StrainData.mat')
for i = 1:343
    m = strains{i};
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

figure
hold on
plot([1:1:343],a(:,3),'k-','LineWidth',2,'color',[215,48,39]/255)
plot([1:1:343],a(:,1),'k-','LineWidth',2,'color','k')
plot([1:1:343],a(:,2),'k-','LineWidth',2,'color',[69,117,180]/255)
set(gca,'FontSize',6,'FontName','Helvetica');
set(gca,'ycolor','k');
set(gca,'xcolor','k');
ylabel('Reaction number','FontSize',8,'FontName','Helvetica','Color','k');
xlabel('Sampled species','FontSize',8,'FontName','Helvetica','Color','k');

set(gcf,'position',[10 50 200 100]);
set(gca,'position',[0.15 0.2 0.65 0.85]);

% generate data for ipath3.0
[~,ID] = ismember(acce_rxn,model_original.rxns);
acce_rxn_kegg = model_original.rxnKEGGID(ID(ID~=0));
[~,ID] = ismember(core_rxn,model_original.rxns);
core_rxn_kegg = model_original.rxnKEGGID(ID(ID~=0));
