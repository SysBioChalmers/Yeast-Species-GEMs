

%Load chemostat data:
fid = fopen('../../data/physiology/Fermentation.tsv','r');
exp_data = textscan(fid,'%s %f32 %f32 %f32 %f32 %f32 %f32 %f32 %s','Delimiter','\t','HeaderLines',1);
species = exp_data{1};
exp_data = [exp_data{2} exp_data{3} exp_data{4} exp_data{5} exp_data{6} exp_data{7} exp_data{8}];
fclose(fid);


currentpath = pwd;

for i = 1:length(exp_data)
    cd(inputpath)
    load([species{i},'.mat'])
    model = reducedModel;
    cd(currentpath)
    mod_data(i) = simulateChemostat(model,exp_data(i,:),1,'C');
end

cd ../modelTests/
% plot the figure
figure
hold on
cols = [215,25,28;253,174,97;171,217,233;44,123,182]/256;
b(1) = plot(exp_data1(:,4),mod_data1(:,4),'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(2,:));
b(2) = plot(exp_data2(:,4),mod_data2(:,4),'s','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(1,:));
b(3) = plot(exp_data3(:,4),mod_data3(:,4),'d','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(3,:));
b(4) = plot(exp_data4(:,4),mod_data4(:,4),'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',cols(4,:));
exp_max = max(exp_data2(:,4));
mod_max = max(mod_data1(:,4));
lim = max(exp_max,mod_max)+0.05;
xlim([0 lim])
ylim([0 lim])
x=0:0.001:lim;
y = x;
plot(x,y,'--','MarkerSize',6,'Color',[64,64,64]/256)
xlabel('Experimental growth rate [1/h]','FontSize',14,'FontName','Helvetica')
ylabel('In silico growth rate [1/h]','FontSize',14,'FontName','Helvetica')
legend(b,'N-limited aerboic','C-limited aerobic','C-limited anaerobic','N-limited anaerobic','Location','northwest')
meanerror = sum(([exp_data1(:,4);exp_data2(:,4);exp_data3(:,4);exp_data4(:,4)]-[mod_data1(:,4);mod_data2(:,4);mod_data3(:,4);mod_data4(:,4)]).^2)/32;
text(0.4,0.1,['meanerror:',num2str(meanerror*100),'%'])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mod_data,solresult] = simulateChemostat(model_origin,exp_data,mode1,mode2)
model = model_origin;
%Relevant positions:
pos(1) = find(strcmp(model.rxns,'r_1714')); %glc
pos(2) = find(strcmp(model.rxns,'r_1992')); %O2
%pos(3) = find(strcmp(model.rxns,'r_1654')); %NH3
pos(3) = find(strcmp(model.rxns,'r_1672')); %co2
pos(4) = find(strcmp(model.rxns,'r_2033')); %pyr
pos(5) = find(strcmp(model.rxns,'r_1761')); %ethanol
pos(6) = find(strcmp(model.rxns,'r_1808')); %gly
pos(7) = find(strcmp(model.rxns,'r_2111')); %growth

%Simulate chemostats:
mod_data = zeros(size(exp_data));
solresult = zeros(length(model.rxns),length(exp_data(:,1)));
if strcmp(mode2,'N')
    model = scaleBioMass(model,'protein',0.289);
    model = scaleBioMass(model,'lipid',0.048);
    model = scaleBioMass(model,'RNA',0.077,'carbohydrate');
end
if mode1 == 2
    model = anaerobicModel(model);
end
for i = 1:length(exp_data(:,1))
    model_test= model;
    %Fix glucose uptake rate and maxmize growth:
    for j = 1:length(exp_data(1,:))-1
        if j <= 2
            model_test = changeRxnBounds(model_test,model_test.rxns(pos(j)),-exp_data(i,j),'l');
        elseif j > 3
           model_test = changeRxnBounds(model_test,model_test.rxns(pos(j)),exp_data(i,j),'l'); % for byproduct production
        end
    end

    model_test = changeObjective(model_test,model_test.rxns(pos(8)),+1);
    sol        = optimizeCbModel(model_test,'max');
    %Store relevant variables:
    mod_data(i,:) = abs(sol.x(pos)');
    solresult(:,i) = sol.x;
end
end
