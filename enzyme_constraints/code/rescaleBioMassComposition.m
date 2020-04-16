 function model = rescaleBioMassComposition(model,phenotype,biomass_phen)
% model = rescaleBioMassComposition(model,Ptot,GAM,scale_comp)
% 
% Benjamin Sanchez. Last update: 2018-10-23
% Ivan Domenzain.   Last update: 2020-04-10

%Get biomass pseudoreaction ID and biomass components pseudoreactions names
cd GECKO/geckomat/limit_proteins
[~,Pbase,Cbase,Rbase,Dbase,Lbase] = sumBioMass(model);
%Get phenotype biomass composition
bio_comp = biomass_phen(:,[1,find(strcmpi(phenotype,biomass_phen.Properties.VariableNames))]);
Ptot     = bio_comp{strcmpi('protein',bio_comp.Type),2};
Ctot     = bio_comp{strcmpi('carbohydrate',bio_comp.Type),2};
Ltot     = bio_comp{strcmpi('lipid',bio_comp.Type),2};
Rtot     = bio_comp{strcmpi('RNA',bio_comp.Type),2};
Dtot     = bio_comp{strcmpi('DNA',bio_comp.Type),2};
GAM      = bio_comp{strcmpi('GAM',bio_comp.Type),2};
%Compute rescaling fractions:
fP = Ptot/Pbase;
fC = Ctot/Cbase;
fL = Ltot/Lbase;
fR = Rtot/Rbase;
fD = Dtot/Dbase;
%Change compositions:
model = rescalePseudoReaction(model,'protein',fP);
model = rescalePseudoReaction(model,'carbohydrate',fC);
model = rescalePseudoReaction(model,'lipid backbone',fL);
model = rescalePseudoReaction(model,'RNA',fR);
model = rescalePseudoReaction(model,'DNA',fD);
%If model contain SLIMER reactions (separate pseudoreactions for
%lipid chains and backbones
model = rescalePseudoReaction(model,'lipid chain',fL);
x= sumBioMass(model);
%Change GAM:
xr_pos = strcmp(model.rxns,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAM = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 & isGAM
        model.S(i,xr_pos) = sign(S_ix)*(GAM);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = rescalePseudoReaction(model,metName,f)
rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
if sum(rxnPos) == 1
    for i = 1:length(model.mets)
        S_ir   = model.S(i,rxnPos);
        isProd = strcmp(model.metNames{i},metName);
        if S_ir ~= 0 && ~isProd
            model.S(i,rxnPos) = f*S_ir;
        end
    end
end
end