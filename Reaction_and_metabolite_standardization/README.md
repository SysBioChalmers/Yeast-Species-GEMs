This file contains R script to analyse model simulation result from 332 yeast species

Three main functions contained in this repo:
1. produce new reaction maxtix
rxnID	coefficient	Metabolite standard name	Metabolite_type	compartment
MNXR96437	1	Ca(2+)	reactant	extracellular
MNXR96437	1	H+	reactant	extracellular
MNXR96437	1	Ca(2+)	product	cytoplasm
MNXR96437	1	H+	product	cytoplasm
MNXR96797	1	chloride	reactant	extracellular
MNXR96797	1	chloride	product	cytoplasm

The above work can be done by 'reaction and metabolite standardization based on mnxID-general.R' if the reaction has MNXID!

2. produce new reaction property
rxnID	rev	GPR	rxn_name_seed	EC	rxnID_kegg	Source;reason
MNXR96437	1	YOL122C	Ca(2+) transport 			bigg:CAt4; gRrules come from TCDB, not included in the DBnewAnnotaiton update (PR #142)

The above work can be done by 'reaction and metabolite standardization based on mnxID-general.R'  if the reaction has MNXID!

3. produce new metabolite property
NewMetName	Charged formula	Charge	compartment	KEGG ID	CHEBI ID	Remark
Ca(2+) [extracellular]	Ca	2	extracellular	C00076	CHEBI:29108	

4. 
However, if the reaction only the ID from kegg database, then we need firstly find the MNXID for each
reaction based on the ID mapping between the kegg ID and MNXID.
The above work can be done by 'reaction and metabolite standardization based on keggID-general.R'  if the reaction has MNXID!



Revised by Hongzhong 2019-8-9