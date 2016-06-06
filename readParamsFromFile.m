%Read parameters from a text file into Matlab workspace.
%Author: Simon Gordonov
%Date: 05/01/2014

function params = readParamsFromFile(paramsFilePath)

    %Read-in the parameter input file (paramsFilePath)
    fid = fopen(paramsFilePath);
    inputParams = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    
    %Remove the comment (starting with '%') and the blank lines in the
    %parameter input file:
    inputParams = inputParams{1};
    locRem1 = cell2mat(cellfun(@(x) isempty(x),inputParams,'UniformOutput',false));
    locRem2 = cellfun(@(x) (regexp(x,'^%')),inputParams,'UniformOutput',false);
    locRem2 = cell2mat(cellfun(@(x) ~isempty(x),locRem2,'UniformOutput',false));
    locRem = locRem1|locRem2;
    inputParams(locRem) = [];
    
    %Set parameters based on input:
    for i = 1:numel(inputParams)
        eval([inputParams{i},';']);
    end
    
    %Save the parameters structure:
    save(fullfile(params.parentFolderForAnalysis,'params.mat'),'params');