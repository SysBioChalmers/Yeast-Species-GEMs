This file contains R script to find the new GPR from panID based on kegg and eggnog web services
if the new KO id occured in kegg database and eggnog, then the annotation from kegg will be used.

Firstly a pipeline is established to conduct the bi-direction blast analyse:
'main_function_BBH.R'
'bi_directional_blast_hit.R'
'standard_bi_direction_blast_analysis.R'

Then, we conduct some initial analysis using:
'analysis of pangenome based reference and non_reference.R'
'analysis of pangenome based on core and variable.R'

Lastly, we get the new GPRs from the annotation from EggNOG and KEGG:
'new_GPR_based_kegg_eggnog.R'


This repo is from Yeast8 project, now used in the new evolution project.
Revised by Hongzhong 2019-8-5


'reciprocal_blast.R' is used to conduct the bi-directional best hit analysis.
Added by Hongzhong 2019-8-7