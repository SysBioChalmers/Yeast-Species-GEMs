function [id,compMat] = compareModelField(models,field)
    % Generates a list of unique field entries and a matrix quantifying the
    % number of appearances of each field entry in each model
    
    % get unique list of field entries
    hasfield = cellfun(@(m) isfield(m,field),models);
    id = cellfun(@(m) m.(field),models(hasfield),'UniformOutput',false);
    id = unique(vertcat(id{:}));
    
    % assemble matrix comparing frequency of each entry in each model
    compMat = zeros(numel(id),numel(models));
    for i = 1:numel(models)
        [~,entryIndex] = ismember(models{i}.(field),id);  % get index of each field entry in the unique id list
        compMat(:,i) = histcounts(entryIndex, 0.5:1:(numel(id)+0.5));  % determine the frequency at which each index appears
    end
end