function [result,gpr] = rxnExistDraft(strains,rxns,type)
% this function is to test whether a rxn exist or not in draft model.
% type here should be 'kegg' or 'metacyc'.
cd ../../../draft_GEM_all_yeast_species/
gpr = cell(1,length(strains)); % create empty cell array for gpr.
for i = 1:length(rxns)
    for j = 1:length(strains)
        if strcmp(type,'metacyc')
        cd 'strain specific model from RAVEN_biocyc_55_110'/
        elseif strcmp(type,'kegg')
        cd 'strain specific model from RAVEN_kegg'/
        else
            error([type,'is not valid. Please enter type: metacyc or kegg'])
        end
        % load draft model
        cd(strains{j})
        fid = fopen('draft_GEM.tsv','r');
        draft = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
        draft = [draft{1} draft{2} draft{3} draft{4} draft{5}];
        fclose(fid);
        cd ../
        % index whether the rxn is in the draft
        [~,idx] = ismember(rxns(i),draft(:,1));
        if ~isempty(idx) && idx ~= 0
            result(i,j) = 1;
            gpr(i,j) = draft(idx,5);
        else
            result(i,j) = 0;
        end
        cd ../
    end
end
cd ../Reconstruction_script/Reconstruction/otherchanges
end
        