%FUNCTION dirSubfolders = readDirSubfolders(inputFolder,expr)
%
% This function finds the subfolders in |inputFolder| and outputs them as a
% cell array of subfolder names in |dirSubfolders|. If |expr| = 'all' then
% all subfolders are outputted, otherwise only subfolder that contain the 
% |expr| substring will be outputted.
%
%Author: Simon Gordonov
%Date: 01/21/14

function dirSubfolders = readDirSubfolders(inputFolder,expr)

    dirSubfolders = dir(inputFolder);
    isDir = [dirSubfolders.isdir];
    dirSubfolders = {dirSubfolders(isDir).name}';
    dirSubfolders(ismember(dirSubfolders,{'.','..'})) = [];
    
    if ~strcmp(expr,'all')
        loc = cell2mat(cellfun(@(x) ~isempty(regexp(x,expr,'once')),...
            dirSubfolders,'UniformOutput',0));
        dirSubfolders = dirSubfolders(loc);
    end
    
