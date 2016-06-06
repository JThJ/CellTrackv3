function [] = cellTrack_ParseStacks_v3(paramsFilePath)

%%%%%%%%
% This code sorts raw image files from live imaging in a folder structure
% and needs to be run only ***ONCE*** as a preperation for analysis
%%%%%%%%
% Copyright MIT 2016
% Developed by Simon Gordonov
% Laboratory for Computational Biology & Biophysics
%%%%%%%%

%% INPUT (PARAMETER FILES AND SCRIPTS)

% What type of file extension are the images
extension = 'tif';

% 1. Read in the parameter file (file paths)
params = readParamsFromFile(paramsFilePath);

% 2. path to all scripts by adding path to the current folder
addpath(genpath(pwd));

% 3. path to new location stacked images is within the current folder with
% the scripts. If it doesn't exist already, a folder "Raw Images" is
% created where corresponding images in each field of view will be
% placed in their own folders.
if ~exist(fullfile(params.parentFolderForAnalysis,'Raw Images'),'dir')
    mkdir(fullfile(params.parentFolderForAnalysis,'Raw Images'));
end

% 4. path to original location images is included in the
% ParseStacks_params.txt file, and can be changed there

% 5. Read in all the raw image file names and store in 'imageFiles':
imageFiles = readDirImages(params.inputImgsPath,extension,0);


%% EXTRACT INFORMATION FILE NAMES

% Extract information in the filenames from each image name
% Make sure to include time and name it (t)

% If you are not familiar with expression notation, leave 'expression' and
% 'Fieldnames' empty. The program 'extractname' will then ask the user to
% 'explain' how the filename is built up, and from that generate an
% expression. After that, the program will loop over the names (e.g.
% time or RowLabel) in 'FieldNames' to create a field in imgInfo with that
% name (e.g. imgInfo.RowLabel)

expression = '^.*_A01_T(\d+)F(\d+)L(\d+)A(\d+)Z(\d+)C(\d+).tif';
FieldNames = {'t'; 'f'; 'l'; 'a'; 'z'; 'c'};

if isempty(expression)
    [expression,FieldNames] = extractname(extension);
end

toks = cellfun(@(x) regexp(x,expression,'tokens'),imageFiles.filenames,...
    'uniformoutput',0);

for i=1:numel(FieldNames)
    imgInfo.(FieldNames{i}) = cellfun(@(x) cat(1,x{1}{i}),toks,'uniformoutput',0);
end

% To make sure that all the files are in the correct order when they are
% read in for subsequent analysis, we need to add leading zeros to the time
% labels. With all other FieldNames we can delete the leading zeros.

for CorrectLeadingZeros=1
    
    for i=1:numel(FieldNames) % loop over all FieldNames
        
        % loop over all the images
        for j=1:numel(toks)
            
            % while the first chr is '0' delete this '0'
            while strcmp(imgInfo.(FieldNames{i}){j,1}(1),'0')
                imgInfo.(FieldNames{i}){j,1} = imgInfo.(FieldNames{i}){j,1}(2:end);
            end
        end
    end
    
    % make a seperate cell field with the new names of the files
    imageFiles.NewFileNames = cell(numel(FieldNames),1);
    
    % determine the length of the latest timepoint
    MaxLength = max(cellfun('length',imgInfo.t));
    
    % loop over all images
    for i=1:numel(toks)
        
        NewFileName = imgInfo.t{i};
        
        % while the string is not long enough, add leading zero
        while length(NewFileName)<MaxLength
            NewFileName = strcat('0',NewFileName);
        end
        
        % save the resulting NewFileName in our list
        imageFiles.NewFileNames{i} = NewFileName;
    end
    
    
end

%% CREATE FOLDER STRUCTURE

% Fill in what FieldNames you would like to use for folder ordering. The
% order of the FieldNames will also govern the mainfolder-subfolder
% structure. Do in this example 'channel1' and 'channel2' will become
% subfolders of the folder 'field1'
FolderStructure = { 'f' , 'field' ;...
    'c' , 'channel'};

% This code creates the appropriate folder structure according to the above
% added FieldNames in FolderStructure
[h,~] = size(FolderStructure);
for i=1:h
    
    if i==1
        uni1 = unique(imgInfo.(FolderStructure{i,1}));
        FolderStructure{i,3} = length(uni1);
        FolderStructure{i,4} = {};
        for j=1:FolderStructure{i,3}
            FolderStructure{i,4}{j} = strcat(FolderStructure{i,2},uni1(j));
            mkdir('Raw Images',char(FolderStructure{i,4}{j}));
        end
    else
        uni2 = unique(imgInfo.(FolderStructure{i,1}));
        FolderStructure{i,3} = length(uni2);
        for j=1:FolderStructure{i,3}
            FolderStructure{i,4}{j} = strcat(FolderStructure{i,2},uni2(j));
            for k=1:length(FolderStructure{i-1,4})
                mkdir(strcat('Raw Images\',char(FolderStructure{i-1,4}{k})),char(FolderStructure{i,4}{j}));
            end
        end
    end
end

%% COPY IMAGES TO FOLDER STRUCTURE

% loop over all the images
for i = 1:numel(toks)
    
    % select the name of the folder (e.g. 'field')
    Folder_loc = strcat(FolderStructure{1,2},imgInfo.(FolderStructure{1,1}){i});
    
    % if necessary add the subfolder path (e.g. 'channel')
    if h==2
        Folder_loc = strcat(Folder_loc,'\',FolderStructure{2,2},imgInfo.(FolderStructure{2,1}){i});
    end
    
    % if necessary add the sub-subfolder path (e.g. 'z-stack')
    if h==3
        Folder_loc = strcat(Folder_loc,'\',FolderStructure{3,2},imgInfo.(FolderStructure{3,1}){i});
    end
    
    % select the course of the original image
    source = imageFiles.filenamesWithPath{i};
    
    % concatenate the target (sub)folder name
    target = strcat('Raw Images\',Folder_loc,'\',imageFiles.NewFileNames{i},'.',extension);
    
    % copy file to new location
    copyfile(source,target);
    
end

end