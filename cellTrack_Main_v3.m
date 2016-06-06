%%%%%%%%
% This code analyses the images to extract cell positions and return them
% in a boolian location array
%%%%%%%%
% Copyright MIT 2016
% Developed by Simon Gordonov
% Laboratory for Computational Biology & Biophysics
%%%%%%%%

function cellTrack_Main_v3(paramsFilePath) %this is cellTrack_Main_params
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT - RUN FOR ANY OF THE BELOW OPTIONS!

    %Add path of the folder that contains all relevant scripts
    addpath(genpath(pwd));
    %addpath(genpath('C:\Users\Simon\Dropbox (MIT)\Lauffenburger Research 2012\Image Processing Scripts\common'))

    %Read in the parameter file 
    params = readParamsFromFile(paramsFilePath);

    %Create output folder path if it doesn't exist
    if ~exist(params.outputPath,'dir')
        mkdir(params.outputPath);
    end
    
    %Save the Main params:
    paramsMain = params;
    save(fullfile(params.parentFolderForAnalysis,'paramsMain.mat'),'paramsMain')
    
    %Read in the folder where frames are stored for all fields:
    fieldFolders = readDirSubfolders(params.inputImgsPath,'all');    

%% Automatically segment the nuclei
    
    %Start parallel pool
    parpool()

    for fieldNum = 1:numel(fieldFolders)
       
        %Get the name of the input folder without path.
        currFieldFolder = fieldFolders{fieldNum};

        %Read in images in current folder, channel 2 should ALWAYS be the
        %nuclei channel for tracking.
        imageFiles.nuc = readDirImages(fullfile(params.inputImgsPath,currFieldFolder,'c2'),'tif',1);
        
        %Check if there is the correct number of time points imaged for
        %given field; if not, continue to next field and move folder
        if numel(imageFiles.nuc.filenames) ~= params.numTimePts
            if ~exist(fullfile(params.parentFolderForAnalysis,'Removed Fields'),'dir')
                mkdir(fullfile(params.parentFolderForAnalysis,'Removed Fields'))
            end
            movefile(fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder),...
                fullfile(params.parentFolderForAnalysis,'Removed Fields',currFieldFolder));
            continue;
        end
        
        %Run the nuclear segmentation (without user- editing). The idea
        %is that the segmentation runs slowly compared to other steps, so it's
        %best if the user is not waiting for the segmentation steps.
        tic
        nucSegData = cell(params.numTimePts,1);
        InucStack = zeros(params.sizeI(1),params.sizeI(2),params.numTimePts);
        parfor currFrameNum = 1:params.numTimePts
            InucStack(:,:,currFrameNum) = mat2gray(imread(imageFiles.nuc.filenamesWithPath{currFrameNum}));
            nucSegData{currFrameNum} = cellTrack_segmentNuclei_v2(InucStack(:,:,currFrameNum),params);
            fprintf('\n%s',['NUCLEI SEGMENTATION FOR FRAME ',num2str(currFrameNum),' OF ',num2str(params.numTimePts),' COMPLETE.']); 
        end
        processedFrames.nucMask = nucSegData;
        toc
        if ~exist(fullfile(params.outputPath,currFieldFolder),'dir')
            mkdir(fullfile(params.outputPath,currFieldFolder));
        end
        
        fprintf('\n%s',['NUCLEI SEGMENTATION FOR FIELD ',num2str(fieldNum),' OF ',num2str(params.numFields),' COMPLETE.']); 
        
        %Save nuclei detections for all fields
        save(fullfile(params.outputPath,currFieldFolder,'processedFrames.mat'),'processedFrames');

    end
    
    %Close parallel pool
    delete(gcp)
    
    fprintf('\n%s\n','NUCLEAR SEGMENTATION COMPLETE.')
    
    %Display nuclei detection for current field (for QC purposes before running detection on all fields):
    InucMaskStack = false(params.sizeI(1),params.sizeI(2),params.numTimePts);
    for currFrameNum = 1:params.numTimePts
        InucMaskStack(:,:,currFrameNum) = imdilate(logical(full(processedFrames.nucMask{currFrameNum})),strel('disk',3,0));
    end
    imseriesshow_multichannel({InucStack,InucMaskStack},'displayRanges',[[0 1];[0 1]])

    %Create new stack based on tracking filtering:
%     InucMaskStackNew = zeros(sizeI(1),sizeI(2),params.numTimePts);
%      for currFrameNum = 1:params.numTimePts
%         tempLocs = vertcat(trackNucCentMtx{currFrameNum,:}); 
%         tempLocs = sub2ind(sizeI,tempLocs(:,2),tempLocs(:,1));
%         tempFrame = zeros(sizeI(1),sizeI(2));
%         tempFrame(tempLocs) = 1;
%         InucMaskStackNew(:,:,currFrameNum) = tempFrame;
%         InucMaskStackNew(:,:,currFrameNum) = imdilate(InucMaskStackNew(:,:,currFrameNum),ones(5,5));
%      end
%         imseriesshow_multichannel({InucStack,InucMaskStackNew},'displayRanges',[[0 1];[0 1]])

% % % %     %TEMP for 122815
% % % %     
% % % %     tNums = cellfun(@(x) regexp(x,'t(\d+)c\d.tif','tokens'),imageFiles.filenames,'uniformoutput',0);
% % % %     tNums = unique(cellfun(@(x) x{1}{1},tNums,'uniformoutput',0));
% % % % 
% % % %     [~,idx] = sort(cell2mat(cellfun(@(x) str2num(x),tNums,'uniformoutput',0)));
% % % %     for i = 2:numel(fieldFolders)
% % % %         currFieldFolder = fieldFolders{i};
% % % %         name = fullfile(params.outputPath,currFieldFolder,'processedFrames.mat');
% % % %         load(name);
% % % %         processedFrames.nucMask = processedFrames.nucMask(idx);
% % % %         save(fullfile(params.outputPath,currFieldFolder,'processedFrames.mat'),'processedFrames');
% % % %     end

%% Setting the cell body segmentation parameters by the user.
% params.cbSeg.numManualsegFrames gives the number of frames
% when user specifies the segmentation threshold. The larger the number,
% the more accurate the segmentation will be in principle over the frames.
% The frames numbers sampled are evenly spaced over the sequence.

%Read in the folder where frames are stored for all fields, but those that
%have frames for all time points after nuclei localization.
fieldFolders = readDirSubfolders(params.outputPath,'all');

if params.boolAutoCBseg == 1

    for fieldNum = 1:numel(fieldFolders)
        
        fprintf('\n%s\n','STARTING SEGMENTATION FOR NEXT FIELD. Please wait...')

        %Get the name of the input folder without path.
        currFieldFolder = fieldFolders{fieldNum};
        
        %Load the stored segmented nuclei masks.
        load(fullfile(params.outputPath,currFieldFolder,'processedFrames.mat'));

        %Read in cell body + nuclei images in current folder, channel 2 should ALWAYS be the
        %nuclei channel for tracking.
        imageFiles.cb = readDirImages(fullfile(params.inputImgsPath,currFieldFolder,'c1'),'tif',0);
        imageFiles.nuc = readDirImages(fullfile(params.inputImgsPath,currFieldFolder,'c2'),'tif',0);

        %Run the cell body segmentation to set thresholds
        frameNumsToSeg = round(linspace(1,params.numTimePts,params.cbSeg.numManualsegFrames));
        if params.cbSeg.numManualsegFrames == 1
            frameNumsToSeg = 1;
        end
        
        %Read in the frames:
        inputFramesCBforSeg = cellfun(@(x) mat2gray(imread(x)),cellstr(vertcat(imageFiles.cb.filenamesWithPath{frameNumsToSeg})),'uniformoutput',0);
        inputFramesNucforSeg = cellfun(@(x) mat2gray(imread(x)),cellstr(vertcat(imageFiles.nuc.filenamesWithPath{frameNumsToSeg})),'uniformoutput',0);

        %Initialize where cell body segmentation thresholds will be saved
        processedFrames.TotsuCB = zeros(params.numTimePts,1);
        processedFrames.TedgeCB = cell(params.numTimePts,1);
        processedFrames.threshSetDescriptors = cell(params.numTimePts,1);

        for currFrameNum = 1:params.cbSeg.numManualsegFrames
            
            Icb = inputFramesCBforSeg{currFrameNum};
            Inuc = inputFramesNucforSeg{currFrameNum};
            InucMask = full(processedFrames.nucMask{frameNumsToSeg(currFrameNum)});
            
             %Run initial segmentation on actin channel, or use parameters
             %from previous segmentation.
            if currFrameNum == 1 
                %Get thresholds for cell body segmentation.
                initialTotsuCB = graythresh(Icb)-0.7*(graythresh(Icb));
                [~,initialTedgeCB] = edge(Icb,'canny');
            else
                initialTotsuCB = processedFrames.TotsuCB(frameNumsToSeg(currFrameNum-1));
                initialTedgeCB = processedFrames.TedgeCB{frameNumsToSeg(currFrameNum-1)};
            end
            
            [processedFrames.TotsuCB(frameNumsToSeg(currFrameNum)),...
                processedFrames.TedgeCB{frameNumsToSeg(currFrameNum)}] = ...
                cellTrack_cb_SegmentThreshSet(Icb,Inuc,InucMask,initialTotsuCB,initialTedgeCB,params);
            
            fprintf('\n%s',['CELL BODY SEGMENTATION THRESHOLD SETTING FOR FRAME ',num2str(frameNumsToSeg(currFrameNum)),' OF ',num2str(params.numTimePts),' COMPLETE.']); 

        end
        
        %Fill in segmentation thresholds for the frames that weren't
        %manually set (by finding parameters for closest frame that was).
%         for currFrameNum = 1:params.numTimePts
%             if processedFrames.TotsuCB(currFrameNum) == 0 %if no thresh value set for this frame
%                 [~,loc] = min(abs(frameNumsToSeg - currFrameNum));
%                 processedFrames.TotsuCB(currFrameNum) = processedFrames.TotsuCB(frameNumsToSeg(loc));
%                 processedFrames.TedgeCB{currFrameNum} = processedFrames.TedgeCB{frameNumsToSeg(loc)};
%             end
%         end
        
        if ~exist(fullfile(params.outputPath,currFieldFolder),'dir')
            mkdir(fullfile(params.outputPath,currFieldFolder));
        end
        
        save(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBThreshSet.mat'),'processedFrames');
        
        close all;
        
        fprintf('\n%s\n',['AUTOMATIC SEGMENTATION THRESHOLD SETTING FOR FIELD ',num2str(fieldNum),' of ',num2str(numel(fieldFolders)),' COMPLETE.']);

    end
end
    
%% Segmentation and object identification (automatrically or manually)

fieldFolders = readDirSubfolders(params.outputPath,'all');

if params.boolAutoCBseg == 1

    for fieldNum = 1:numel(fieldFolders)

        %Get the name of the input folder without path.
        currFieldFolder = fieldFolders{fieldNum};

        %Load the stored segmented nuclei masks.
%         load(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBThreshSet.mat'));
        load(fullfile(params.outputPath,currFieldFolder,'processedFrames.mat'));

%         %TEMP!!!!!!!!!!!!
%         processedFrames.TedgeCB(2:end) = cell(numel(processedFrames.TedgeCB)-1,1); 
%         processedFrames.TotsuCB(2:end) = 0;
        
        %Get number of frames
        numFrames = numel(processedFrames.nucMask);
        
        %Read in images in current folder, channel 2 should ALWAYS be the
        %nuclei channel for tracking.
        imageFiles.nuc = readDirImages(fullfile(params.inputImgsPath,currFieldFolder,'c2'),'tif',0);
        imageFiles.cb = readDirImages(fullfile(params.inputImgsPath,currFieldFolder,'c1'),'tif',0);
        
        inputFramesNuc = cellfun(@(x) imread(x),... 
                        imageFiles.nuc.filenamesWithPath,'UniformOutput',0);
        inputFramesCB = cellfun(@(x) imread(x),... 
                        imageFiles.cb.filenamesWithPath,'UniformOutput',0);
        
        %%%%%%
        processedFrames.refTocurrFrameSegScoreMtx = cell(numFrames,1);
        processedFrames.TotsuVec = cell(numFrames,1);
        processedFrames.TedgeVec = cell(numFrames,1);
        processedFrames.cbMaskFieldLocal = cell(numFrames,1);
        processedFrames.cbMaskFieldLocalSplit = cell(numFrames,1);
        processedFrames.nucMaskFieldLocal = cell(numFrames,1); 
        %%%%%%
        
        %Store the cell body images in stack array
        Icb = mat2gray(imread(imageFiles.cb.filenamesWithPath{1}));
        processedFrames.frameSize = size(Icb);
        IcbStack = zeros(size(Icb,1),size(Icb,2),params.numTimePts);
        
        for currFrameNum = 1:numFrames

            %Read in and create an array of current folder frame sequence.
            Icb = mat2gray(inputFramesCB{currFrameNum});
            IcbStack(:,:,currFrameNum) = Icb;
            
            %Set the nuclear segmentation into processed array done previously
            %(automatically).
            processedFrames.maskNuc{currFrameNum} = logical(full(processedFrames.nucMask{currFrameNum}));
            processedFrames.maskNucShrink{currFrameNum} = bwmorph(processedFrames.maskNuc{currFrameNum},'shrink',Inf);

            %Run cell body segmentation automatically using the cell body
            %threshold parameters set by user (and interpolated into
            %intermediate frames) in preceding steps.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            
%             prevProcessedFrames = load(fullfile('Y:\homedirs\gordonov\Mun_Dox_231LAH2B_2DMigr_120414\Output\Tracking Data',...
%                 currFieldFolder,'processedFramesWithCBThreshSet.mat'));
%             prevProcessedFrames = prevProcessedFrames.processedFrames;
%             
%             processedFrames.TotsuCB(currFrameNum) = prevProcessedFrames.TotsuCB(currFrameNum);
%             processedFrames.TedgeCB(currFrameNum) = prevProcessedFrames.TedgeCB(currFrameNum);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             For 010815 experiment ONLY, load the thresholds set before (not v2):
            prevProcessedFrames = load(fullfile('Y:\homedirs\gordonov\Mun_DoxAndSB_231LAH2B_2DMigr_010815\Output\Tracking Data',...
                currFieldFolder,'processedFramesWithCBThreshSet.mat'));
            prevProcessedFrames = prevProcessedFrames.processedFrames;
            
            processedFrames.TotsuCB(currFrameNum) = prevProcessedFrames.TotsuCB(currFrameNum);
            processedFrames.TedgeCB(currFrameNum) = prevProcessedFrames.TedgeCB(currFrameNum);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %If the Otsu and Canny thresholds for current frame were
            %manually set, just run segmentation using the set thresholds.
            %Otherwise, find thresholds that incrementally preserve
            %descriptor values from nearby frame.
            if processedFrames.TotsuCB(currFrameNum) == 0
                manualSetBool = 0;
            else
                manualSetBool = 1;
            end
           
            processedFrames = cellTrack_cb_SegmentAutoFromThreshSet(processedFrames,Icb,params,manualSetBool,currFrameNum); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %THIS IS THE NEW ADDED PART 03/16/16 - after global segmentation,
            %perform local segmentation per object to "refine" the masks for
            %the cell bodies and nuclei.
            IcbNoTophat = inputFramesCB{currFrameNum};
            InucNoTophat = inputFramesNuc{currFrameNum};
            [processedFrames.cbMaskFieldLocal{currFrameNum},...
                processedFrames.nucMaskFieldLocal{currFrameNum}]...
                = cellTrack_segmentCBandNucmaskFieldLocal_v2(processedFrames,...
                IcbNoTophat,InucNoTophat,params,currFrameNum);
            
            %Make all masks sparse for saving 
            processedFrames.cbMask{currFrameNum} = sparse(processedFrames.cbMask{currFrameNum});
            processedFrames.nucMask{currFrameNum} = sparse(processedFrames.nucMask{currFrameNum});
            processedFrames.nucMaskFieldLocal{currFrameNum} = sparse(processedFrames.nucMaskFieldLocal{currFrameNum});
            processedFrames.cbMaskFieldLocal{currFrameNum} = sparse(processedFrames.cbMaskFieldLocal{currFrameNum});
            processedFrames.maskNucShrink{currFrameNum} = sparse(processedFrames.maskNucShrink{currFrameNum});

            %Notification of segmentation of current frame    
            fprintf('\n%s',['RUNNING SEGMENTATION MODULE FOR FRAME ',...
                num2str(currFrameNum),'/',num2str(numFrames),' of field ',num2str(fieldNum),'/',num2str(numel(fieldFolders)),' COMPLETE.'])
                        
        end
        
        processedFrames = rmfield(processedFrames,'maskNuc');
        
        %Save the processed frames temporary data (nuclear + cell body
        %segmentations)
        save(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBAutoSeg.mat'),'processedFrames');
        
    end

else

    %Desginate current field folder
    currFieldFolder = fieldFolders{fieldNum};

    %Load the stored segmented nuclei masks.
    load(fullfile(params.outputPath,currFieldFolder,'processedFrames.mat'));

    %Get number of frames
    numFrames = numel(processedFrames.nucMask);

    %Define the temporary cell array to store processing for frames in this
    %field of view.
    processedFramesTemp = cell(numFrames,1);

    %Read in images in current folder, channel 2 should ALWAYS be the
    %nuclei channel for tracking.
    imageFiles.nuc = readDirImages(fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder,'c2'),'tif',0);
    imageFiles.cb = readDirImages(fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder,'c1'),'tif',0);

    %For the first frame, user sets the cell body threshold for
    %segmentation, then it's done automatically. It's usually better to
    %make the threshold slightly lower so that we'd rather get cells
    %touching than single cells have improparly segmented regions due to
    %decrease in intensity.

        %Read in and create an array of current folder frame sequence.
        inputFrameCB =  imread(imageFiles.cb.filenamesWithPath{currFrameNum});
        Icb = im2double(inputFrameCB);
        inputFrameNuc = imread(imageFiles.nuc.filenamesWithPath{currFrameNum});
        Inuc = im2double(inputFrameNuc);
        processedFramesTemp{currFrameNum}.frameSize = size(Inuc);

        %Run initial segmentation on nuclei and actin channel.
        if currFrameNum == 1 
            %Get thresholds for cell body segmentation.
            initialTotsuCB = graythresh(Icb)-0.4*(graythresh(Icb));
            [~,initialTedgeCB] = edge(Icb,'canny');
            fprintf('\n%s\n',['ANALYSIS FOR ',currFieldFolder,':']);
        else
            initialTotsuCB = processedFramesTemp{currFrameNum-1}.cbTotsuFinal;
            initialTedgeCB = processedFramesTemp{currFrameNum-1}.cbTedgeFinal;
        end

        %Run segmentation module.
        fprintf('\n%s%d%s%d\n','STARTING SEGMENTATION MODULE FOR FRAME ',...
            currFrameNum,' OF ',numFrames)

        %Run interactive segmentation on the nuclei image.  Note that
        %nuclei segmentation now runs without user-updated thresholds, but
        %only allows to add, remove, and split regions (see
        %cellTrack_nuc_Segment.m).
        InucMask = full(processedFrames.nucMask{currFrameNum});
        processedFramesTemp{currFrameNum} = ...
            cellTrack_nuc_Segment(processedFramesTemp{currFrameNum},Inuc,InucMask,Icb,params);

        %Run interactive segmentation on the cell body image. Note that
        %this segmentation process automatically labels for removal (in
        %red) cells that touch each other, as detected by nuclei
        %segmentation, if nuclei channel is measured.
        if currFrameNum == 1 
            prevProcessedFrameTemp = [];
        else
            prevProcessedFrameTemp = processedFramesTemp{currFrameNum-1};
        end
        processedFramesTemp{currFrameNum} = ...
            cellTrack_cb_Segment(processedFramesTemp{currFrameNum},Inuc,Icb,params,...
            initialTotsuCB,initialTedgeCB,prevProcessedFrameTemp);  
        
     %Delete unneeded fields and compress data
        processedFramesTemp{currFrameNum} = rmfield(processedFramesTemp{currFrameNum},{'maskNuc','maskNucShrink'});
        processedFramesTemp{currFrameNum}.cbMask = sparse(processedFramesTemp{currFrameNum}.cbMask);
        processedFramesTemp{currFrameNum}.nucMask = sparse(processedFramesTemp{currFrameNum}.nucMask);
        processedFramesTemp{currFrameNum}.modifiedRegionsCB = sparse(processedFramesTemp{currFrameNum}.modifiedRegionsCB);

        fprintf('\n%s',['MANUAL SEGMENTATION PROCESSING FOR FIELD',currFieldFolder,'of',num2str(numel(fieldFolders)),' COMPLETE.']);

    %Save the processed frames temporary data (nuclear + cell body
    %segmentations)
    save(fullfile(params.outputPath,currFieldFolder,'processedFramesTemp.mat'),'processedFramesTemp');
   
end

%% Visualize segmentation (optional):

fieldNum=1;

%Get the name of the input folder without path.
currFieldFolder = fieldFolders{fieldNum};

%Load the stored segmented nuclei masks.
load(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBAutoSeg.mat'));

%Get number of frames
numFrames = numel(processedFrames.nucMask);

%Read in images in current folder, channel 2 should ALWAYS be the
%nuclei channel for tracking.
imageFiles.nuc = readDirImages(fullfile(params.inputImgsPath,currFieldFolder,'c2'),'tif',0);
imageFiles.cb = readDirImages(fullfile(params.inputImgsPath,currFieldFolder,'c1'),'tif',0);

%Store the cell body images in stack array
Icb = mat2gray(imread(imageFiles.cb.filenamesWithPath{1}));
processedFrames.frameSize = size(Icb);
IcbStack = zeros(size(Icb,1),size(Icb,2),params.numTimePts);

for currFrameNum = 1:numFrames
    %Read in and create an array of current folder frame sequence.
    Icb = mat2gray(imread(imageFiles.cb.filenamesWithPath{currFrameNum}));
    IcbStack(:,:,currFrameNum) = Icb;
end

IcbMaskStack = false(processedFrames.frameSize(1),processedFrames.frameSize(2),params.numTimePts);
InucMaskStack = false(processedFrames.frameSize(1),processedFrames.frameSize(2),params.numTimePts);
for currFrameNum = 1:params.numTimePts
    IcbMaskStack(:,:,currFrameNum) = logical(full(processedFrames.cbMask{currFrameNum}));
    InucMaskStack(:,:,currFrameNum) = logical(full(processedFrames.nucMaskFieldLocal{currFrameNum}));
end
imseriesmaskshow(IcbStack,{IcbMaskStack,InucMaskStack},'displayRange',[0 1],'maskAlphas',[0.2,0.2],'maskColors',[0 1 0;1 0 0])



%% Perform automated tracking and saving of tracks movies prior to editing tracks.

badfields = {};
fieldFolders = readDirSubfolders(params.outputPath,'all');

for fieldNum = 1:numel(fieldFolders)

    %Get the name of the input folder without path.
        currFieldFolder = fieldFolders{fieldNum};

    %Load temporary processed frames segmentation data
        load(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBAutoSeg.mat'));
        
    %Remove nuclei centers that don't occur within cells after local
    %segmentation
        for i = 1:numel(processedFrames.nucMask)
%             processedFrames.cbMask{i} = full(processedFrames.cbMask{i});
%             temp = logical(processedFrames.cbMask{i}) & ~imclearborder(logical(processedFrames.cbMask{i}));
%             processedFrames.cbMask{i}(temp == 1) = 3;
            processedFrames.nucMask{i} = full(processedFrames.nucMask{i}).*logical(full(processedFrames.cbMaskFieldLocal{i}));
        end
        
    %Read in cell body and nuclei images in current folder, channel 2
    %should ALWAYS be the nuclei channel for tracking.
        imageFiles.nuc = readDirImages(fullfile(params.inputImgsPath,...
            currFieldFolder,'c2'),'tif',0);
        imageFiles.cb = readDirImages(fullfile(params.inputImgsPath,...
            currFieldFolder,'c1'),'tif',0);
    
        inputFramesNuc = cellfun(@(x) imread(x),... 
                        imageFiles.nuc.filenamesWithPath,'UniformOutput',0);
        inputFramesCB = cellfun(@(x) imread(x),... 
                        imageFiles.cb.filenamesWithPath,'UniformOutput',0);
        
    %%%%%%%%%%%%%% Label and track individual cell tracks (objects) %%%%%%%%%%%%%
        params.maxDistTrack = 70;  %Set initial value 

        trackData = cellTrack_LabelAndTrack_withNuc_v2(processedFrames,...
            params,inputFramesCB,inputFramesNuc);
        if isempty(trackData)
            disp('Error in tracking current field')
            badfields = [badfields;currFieldFolder];
            xlswrite(fullfile(params.parentFolderForAnalysis,'Bad fields.xlsx'),badfields)
            continue;
        end

    %%%%%%%%%%%%%%% Output movie pre-editing of tracks (this is useful for
    %%%%%%%%%%%%%%% deleting/combining tracks in the GUI %%%%%%%%%%%%%
                    
    %Get number of frames
       numFrames = numel(inputFramesCB);
                    
    %Make the images green if one channel, and if 2 make nuclei
    %magenta and cell body green.
       inputFramesRGB = cell(numFrames,1);
       for i = 1:numFrames
           inputFramesRGB{i} = imfuse(inputFramesCB{i},inputFramesNuc{i});
       end
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %        %Combine the nuclei and actin channels into and RGB image.             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %        for i = 1:numFrames
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %            inputFramesRGB{i} = cat(3,zeros(processedFramesTemp{1}.frameSize),...
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                inputFramesCB{i},inputFramesNuc{i});
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Icb = inputFramesCB{i};
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Icb = imadjust(Icb);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Icb(logical(processedFramesTemp{i}.cbMask) == 0) = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Icb = imcomplement(Icb);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Inuc = inputFramesCB{i};
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Inuc = imadjust(Inuc);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Inuc(logical(processedFramesTemp{i}.nucMask) == 0) = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Inuc = imcomplement(Inuc);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Ibg = (logical(processedFramesTemp{i}.cbMask == 0)).*255;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Icb(Ibg == 255) = 255;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             Inuc(Ibg == 255) = 255;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             inputFramesRGB{i} = cat(3,Ibg,Icb,Inuc);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             inputFrames = inputFramesRGB;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             imshow(inputFramesRGB{i},[]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %        end             

    %Create AVI file of image frames + tracks       
    if strcmp(params.microscope,'NSD')
       timeStampsStr = xlsread(fullfile(params.inputImgsPath,...
           currFieldFolder,['TimeStamps_',currFieldFolder,'.xlsx']));
       timeStampsStr = timeStampsStr + (params.numHoursStimElapse/24);
       timeStampsStr = cellfun(@(x) strsplit(x,' '),cellstr(datestr(timeStampsStr+693960)),'uniformoutput',0);
       timeStampsStr = vertcat(timeStampsStr{:});
       timeStampsStr = timeStampsStr(:,2);
    elseif strcmp(params.microscope,'Cellomics')
       timeStampsStr = textscan(num2str(1:params.numTimePts),'%s');
       timeStampsStr = timeStampsStr{:};
    elseif strcmp(params.microscope,'Incucyte')
        load(fullfile(params.parentFolderForAnalysis,'timeStamps_SecondsPostTreat.mat'));
        %Convert to time string: (h:m:s):
        timeStampsStr = cell(params.numTimePts,1);
        for i = 1:params.numTimePts
            timeCurr = timeStamps_SecondsPostTreat{1}(i);
            hrs = floor(timeCurr/3600);
            mins = floor(((timeCurr/3600)-hrs)/60);
            secs = timeCurr - (hrs*3600) - (mins*60);
            timeStampsStr{i} = cellstr(datestr(datenum([0 0 0 hrs mins secs]),'HH:MM:SS'));
        end
    end
    
   %Set frame rate and quality of movies to output:
   frameRate = 10;
   movQuality = 100;
   pxScale = params.sizePix;
   
   tracksMovieSaveFolder = fullfile(fileparts(params.outputPath),...
     'time lapse movies with tracks pre-edit',[currFieldFolder,'_Lineage.mp4']);
   if ~exist(fileparts(tracksMovieSaveFolder),'dir')
       mkdir(fileparts(tracksMovieSaveFolder));
   end
   cellTrack_createLineageMovie(inputFramesRGB,processedFrames,trackData,...
       tracksMovieSaveFolder,frameRate,movQuality,timeStampsStr,pxScale);
   
   tracksMovieSaveFolder = fullfile(fileparts(params.outputPath),...
     'time lapse movies with tracks pre-edit',[currFieldFolder,'_SegAndTracks.mp4']);
   if ~exist(fileparts(tracksMovieSaveFolder),'dir')
       mkdir(fileparts(tracksMovieSaveFolder));
   end
   cellTrack_createSegAndTrackMovie(inputFramesRGB,processedFrames,trackData,...
       tracksMovieSaveFolder,frameRate,movQuality,timeStampsStr,pxScale);
   
    %Save the tracking data for current field
    trackData.cbMask = cellfun(@(x) sparse(x),trackData.cbMask,'uniformoutput',0);
    trackData.nucCentMask = cellfun(@(x) sparse(x),trackData.nucCentMask,'uniformoutput',0);
    trackData.nucCentMaskAll = cellfun(@(x) sparse(x),trackData.nucCentMaskAll,'uniformoutput',0);
    trackData.nucMask = cellfun(@(x) sparse(x),trackData.nucMask,'uniformoutput',0);
    trackData.divStatus= sparse(trackData.divStatus);
    trackData.daughtCell = sparse(trackData.daughtCell);
    trackData.parentCell = sparse(trackData.parentCell);
    save(fullfile(params.outputPath,currFieldFolder,'trackData.mat'),'trackData');

    close all;
    
end

%%%%%%%% FROM HERE GO TO:
%cellTrack_analyzeTrackData_v2.m

%% Perform below automated tracking if ONLY on nuclei (NO CELL BODY
% SEGMENTATION). In case of this we analyze migration metrics right away
% from the tracks (cellTrack_Analyze_onlyNuc.m)

badfields = {};
fieldFolders = readDirSubfolders(params.outputPath,'all');

for fieldNum = 1:numel(fieldFolders)

    %Get the name of the input folder without path.
        currFieldFolder = fieldFolders{fieldNum};

    %Load temporary processed frames segmentation data
        load(fullfile(params.outputPath,currFieldFolder,'processedFrames.mat'));
        
    %Make matrices full from sparse
        processedFrames.nucMask = cellfun(@(x) full(x), processedFrames.nucMask,'uniformoutput',0);

    %Set frame size
        processedFrames.frameSize = size(processedFrames.nucMask{1});  
        
    %%%%%%%%%%%%%% Label and track individual cell tracks (objects) %%%%%%%%%%%%%
        params.maxDistTrack = 70; 
        while 1 %Excessive combinatorics error may be encountered, so keep decreasing maxDistTrack till no error occurs
            try
               objects = cellTrack_LabelAndTrack_onlyNuc(processedFrames,params);
               break;
            catch
               params.maxDistTrack = params.maxDistTrack - 1;
            end
            
            if params.maxDistTrack < 10
                badfields = [badfields;currFieldFolder];
                xlswrite(fullfile(params.parentFolderForAnalysis,'Bad fields.xlsx'),badfields)
                break;
            end
        end
        
        if params.maxDistTrack < 10
            continue;
        end
        
        %Save the objects
        save(fullfile(params.outputPath,currFieldFolder,'objects.mat'),'objects');
        
        fprintf('\n%s',['Tracking cell in field ',num2str(fieldNum),' of ',num2str(numel(fieldFolders)),' complete.']);
        
        processedFramesForAnalysis.params = params;
        processedFramesForAnalysis.params.inputImgsPath = fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder);
        processedFramesForAnalysis.objects = objects;
        processedFramesForAnalysis.frameSize = processedFrames.frameSize;
        save(fullfile(params.outputPath,currFieldFolder,'processedFramesForAnalysis.mat'),'processedFramesForAnalysis');  
        
        fprintf('\n%s',['SAVING OBJECT AND PROCESSED DATA FOR FIELD ',num2str(fieldNum),' of ',num2str(numel(fieldFolders)),' COMPLETE.']);
         
    %Read in cell body and nuclei images in current folder, channel 2
    %should ALWAYS be the nuclei channel for tracking.
        imageFiles.nuc = readDirImages(fullfile(params.parentFolderForAnalysis,...
            'Raw Images',currFieldFolder,'c2'),'tif',0);
        imageFiles.cb = readDirImages(fullfile(params.parentFolderForAnalysis,...
            'Raw Images',currFieldFolder,'c1'),'tif',0);
    
        inputFramesNuc = cellfun(@(x) mat2gray(imread(x)),... 
                        imageFiles.nuc.filenamesWithPath,'UniformOutput',0);
        inputFramesCB = cellfun(@(x) mat2gray(imread(x)),... 
                        imageFiles.cb.filenamesWithPath,'UniformOutput',0);               
     %Get number of frames
        numFrames = numel(inputFramesCB);
                    
     %Make the images green if one channel, and if 2 make nuclei
     %magenta and cell body green.
         inputFramesRGB = cell(numFrames,1);
         for i = 1:numFrames
            inputFramesRGB{i} = imfuse(inputFramesCB{i},inputFramesNuc{i});
         end

    %Create AVI file of image frames + tracks       
    if strcmp(params.microscope,'NSD')
       timeStampsStr = xlsread(fullfile(params.parentFolderForAnalysis,...
           'Raw Images',currFieldFolder,['TimeStamps_',currFieldFolder,'.xlsx']));
       timeStampsStr = cellfun(@(x) strsplit(x,' '),cellstr(datestr(timeStampsStr+693960)),'uniformoutput',0);
       timeStampsStr = vertcat(timeStampsStr{:});
       timeStampsStr = timeStampsStr(:,2);
    elseif strcmp(params.microscope,'Cellomics')
       timeStampsStr = textscan(num2str(1:params.numTimePts),'%s');
       timeStampsStr = timeStampsStr{:};
    elseif strcmp(params.microscope,'Incucyte')
        load(fullfile(params.parentFolderForAnalysis,'timeStamps_SecondsPostTreat.mat'));
        %Convert to time string: (h:m:s):
        timeStampsStr = cell(params.numTimePts,1);
        for i = 1:params.numTimePts
            timeCurr = timeStamps_SecondsPostTreat{1}(i);
            hrs = floor(timeCurr/3600);
            mins = floor(((timeCurr/3600)-hrs)/60);
            secs = timeCurr - (hrs*3600) - (mins*60);
            timeStampsStr{i} = cellstr(datestr(datenum([0 0 0 hrs mins secs]),'HH:MM:SS'));
        end
    end
    
   %Set frame rate and quality of movies to output:
   frameRate = 10;
   movQuality = 100;
   pxScale = 1.3;
   tracksMovieSaveFolder = fullfile(fileparts(params.outputPath),...
     'time lapse movies with tracks pre-edit',[currFieldFolder,'_nucTracks.mp4']);
   if ~exist(fileparts(tracksMovieSaveFolder),'dir')
       mkdir(fileparts(tracksMovieSaveFolder));
   end

   createAVIfromImgSeqWithTracks_onlyNuc(inputFramesRGB,objects,tracksMovieSaveFolder,frameRate,movQuality,timeStampsStr,pxScale);

   close all;
end

%%

%ANALYZE MIGRATION BASED ON NUCLEI ONLY

%Analyze migration metrics for tracked nuclei (this is when we DON'T care
%about touching or dividing cells, but just all tracked nuclei in the
%population).
paramsAnalyze = readParamsFromFile(fullfile(params.parentFolderForAnalysis,'cellTrack_Analyze_params.txt'));
paramsAnalyze.condLabs = {'treat'};
paramsAnalyze.condRegExpString = '^(\S+)_[0-9]'; 
%Set the time scale of number of steps to compute displacement,
%total distance travelled, speeds, etc. of the cell centroids.

maxNumNaNsAllowed = 85; %remove a cell if it wasn't tracked for maxNumNaNs or more time points (NaNs in data)?
plot95CIbool = 0; % plot (1) or dont (0) 95% CI for the means in the plots
tScales = [6,12,24]; %if image at 10-minute intervals then tScales = [6,12,24] is 1hr, 2hrs, 4hrs tScales
% condOrder = [7,10,11,9,8,4,3,6,2,1,5]; %set order of conditions for plotting
% condGroups = [1,1,1,2,2,3,3,3,4,4,4]; %set the grouping of conditions for colors
% condCtrl = [1,1,1,4,4,6,6,6,9,9,9];%note this order of controls for the treated conditions AFTER re-ordering of conditions by condOrder.
condOrder = [7,8,9,10,11,12,4,5,6,1,2,3]; %set order of conditions for plotting
condGroups = [1,1,1,2,2,2,3,3,3,4,4,4]; %set the grouping of conditions for colors - note this order is BEFORE re-ordering of conditions by condOrder
condCtrl = [1,1,1,4,4,4,7,7,7,10,10,10];%note this order of controls for the treated conditions AFTER re-ordering of conditions by condOrder.
for i = 1:numel(tScales)
    tScale = tScales(i); %number of frame transitions (steps) e.g.  consider sequence o-o-o; is has tScale of 2 not 3, where o is time point and - is steps (tScale) between time points.
    cellTrack_Analyze_onlyNuc(paramsAnalyze,tScale,maxNumNaNsAllowed,plot95CIbool,params.sizePx,condOrder,condGroups,condTreat,condCtrl);
end

%%
    %MANUSCRIPT SUPPLEMENTAL MATERIAL VIDEOS

    folds = [17,58,68];
    
    for foldNum = 1:numel(folds)
    currFieldFolder = fieldFolders{folds(foldNum)};
    load(fullfile(params.outputPath,currFieldFolder,'objects.mat'),'objects');
    load(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBAutoSeg.mat'));
    tFrame = 48;
    scalePix = 1.3; %microns per pixel
   
    imageFiles.nuc = readDirImages(fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder,'c2'),'tif',0);
        imageFiles.cb = readDirImages(fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder,'c1'),'tif',0);
    
    %%%%%%%%%%%%%%% Output movie pre-editing of tracks (this is useful for
    %%%%%%%%%%%%%%% deleting/combining tracks in the GUI %%%%%%%%%%%%%
       inputFramesNuc = cellfun(@(x) mat2gray(imread(x)),... 
                        imageFiles.nuc.filenamesWithPath,'UniformOutput',0);
       inputFramesCB = cellfun(@(x) mat2gray(imread(x)),... 
                        imageFiles.cb.filenamesWithPath,'UniformOutput',0);
                    
       %Get number of frames
       numFrames = numel(inputFramesCB);
                    
     %Make the images green if one channel, and if 2 make nuclei
     %magenta and cell body green.
         inputFramesRGB = cell(numFrames,1);
         for i = 1:numFrames
            inputFramesRGB{i} = imfuse(inputFramesCB{i},inputFramesNuc{i});
         end

    %Create AVI file of image frames + tracks       
    if strcmp(params.microscope,'NSD')
       timeStampsStr = xlsread(fullfile(params.parentFolderForAnalysis,...
           'Raw Images',currFieldFolder,['TimeStamps_',currFieldFolder,'.xlsx']));
       timeStampsStr = cellfun(@(x) strsplit(x,' '),cellstr(datestr(timeStampsStr+693960)),'uniformoutput',0);
       timeStampsStr = vertcat(timeStampsStr{:});
       timeStampsStr = timeStampsStr(:,2);
    elseif strcmp(params.microscope,'Cellomics')
       timeStampsStr = textscan(num2str(1:params.numTimePts),'%s');
       timeStampsStr = timeStampsStr{:};
    elseif strcmp(params.microscope,'Incucyte')
        load(fullfile(params.parentFolderForAnalysis,'timeStamps_SecondsPostTreat.mat'));
        %Convert to time string: (h:m:s):
        timeStampsStr = cell(params.numTimePts,1);
        for i = 1:params.numTimePts
            timeCurr = timeStamps_SecondsPostTreat{1}(i);
            hrs = floor(timeCurr/3600);
            mins = floor(((timeCurr/3600)-hrs)/60);
            secs = timeCurr - (hrs*3600) - (mins*60);
            timeStampsStr{i} = cellstr(datestr(datenum([0 0 0 hrs mins secs]),'HH:MM:SS'));
        end
    end
    
   %Set frame rate and quality of movies to output:
   frameRate = 10;
   movQuality = 100;
   pxScale = 1.3;
   tracksMovieSaveFolder = fullfile(fileparts(params.outputPath),...
     'time lapse movies with tracks pre-edit',[currFieldFolder,'_MANUSCRIPT_VIDEO']);
   if ~exist(fileparts(tracksMovieSaveFolder),'dir')
       mkdir(fileparts(tracksMovieSaveFolder));
   end
   
   createAVIfromImgSeqWithTracks_withNuc(inputFramesRGB,objects,tracksMovieSaveFolder,frameRate,movQuality,timeStampsStr,pxScale);

   close all;
    end
    
%%
    %MANUSCRIPT FIGURE 2C
    %The script below enables us to plot, for a given field and time frame
    %number, the nuclei tracks up to the selected time frame and segmented
    %cells outlined. The user can select the colormap to use for the numcei
    %tracks.
    currFieldFolder = fieldFolders{14};
    load(fullfile(params.outputPath,currFieldFolder,'objects.mat'),'objects');
    load(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBAutoSeg.mat'));
    tFrame = 48;
    scalePix = 1.3; %microns per pixel
    
    %Read in images in current folder, channel 2 should ALWAYS be the
    %nuclei channel for tracking.
        imageFiles.nuc = readDirImages(fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder,'c2'),'tif',0);
        imageFiles.cb = readDirImages(fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder,'c1'),'tif',0);
        
        cbMask = processedFrames.cbMask{tFrame};

    %%%%%%%%%%%%%%% Output movie pre-editing of tracks (this is useful for
    %%%%%%%%%%%%%%% deleting/combining tracks in the GUI %%%%%%%%%%%%%
       inputFramesNuc = mat2gray(imread(imageFiles.nuc.filenamesWithPath{tFrame}));
       inputFramesCB = mat2gray(imread(imageFiles.cb.filenamesWithPath{tFrame}));

     %Make the images green if one channel, and if 2 make nuclei
     %magenta and cell body green.
%          inputFramesRGB = imfuse(inputFramesCB,inputFramesNuc);
        inputFramesRGB = imfuse(imadjust(inputFramesCB).*full(cbMask),inputFramesNuc.*full(cbMask));
    
    cmapTracks = gray(256);%labmap([pi/3 pi/3],[60 80],[0 60]);%parula(tFrame);
    cmapTracks = cmapTracks(round(linspace(50,size(cmapTracks,1)-100,tFrame-1)),:);
    showTimeFrameSegCellsAndNucleiTracks(inputFramesRGB,objects,cmapTracks,tFrame,scalePix)
    export_fig('C:\Users\Simon\Dropbox (MIT)\Lauffenburger Research 2012\Bayesian, HMM, Regression, PLS Image Data Pilot Project\MBoC Figures\Fig2C.pdf','-painters','-transparent','-append') 

%    fprintf('\n%s',['OBJECTS TRACKING/LABELING AND PRE-EDIT MOVIE SAVING FOR FIELD ',num2str(fieldNum),' of ',num2str(numel(fieldFolders)),' COMPLETE.']);
%%

        
%%EDIT OUT TRACK QC/SELECTION FOR NOW:

% %%  Perform editing of tracks
%         
%     %%%%%%%%%%%%%%% Enable user to delete, combine, and split tracks %%%%%%%%%%%%%
%         cellTrack_EditTracks_withNuc(processedFramesTemp,objects);  
%         
% %% Output tracks movie after editing
% 
%        editedObjects = data.objects; %note data struct come from assigning data to workspace within cellTrack_EditTracks_withNuc subfunctions
%    
%        %%%%%%%%%%%%%%%% Create movie of cells with tracks overlaid
%        tracksMovieSaveFolder = fullfile(fileparts(params.outputPath),...
%          'time lapse movies with tracks post-edit',currFieldFolder);
%        if ~exist(fileparts(tracksMovieSaveFolder),'dir')
%            mkdir(fileparts(tracksMovieSaveFolder));
%        end
%        createAVIfromImgSeqWithTracks_withNuc(inputFramesRGB,editedObjects,tracksMovieSaveFolder...
%            ,frameRate,movQuality,timeStampsStr);
%% 
% Place all data needed for subsequent analysis and graphing of motility/
% morphology property dynamics in frames data and output. Then move on the
% the cellTrack_Analyze_with_Nuc script.
% For the PCA-HMM work, after running code below move on to the
% RegressModel_main script, where we consolidate data, select cells that
% have minimum number of consecutive tracked frames without touching other
% cells or the boundaries, manual segmentation QC GUI, computing
% descriptors, doing PCA after normalization/standardization, running HMM,
% and doing all HMM-based calculations.

%Loop through all the fields for which cells were tracked
fieldFolders = readDirSubfolders(params.outputPath,'all');

for fieldNum = 1:numel(fieldFolders)
    
        %Get the name of the input folder without path.
        currFieldFolder = fieldFolders{fieldNum};
       
        %Load the data generated from various steps in segmentation and
        %tracking:
        load(fullfile(params.outputPath,currFieldFolder,'processedFramesWithCBAutoSeg.mat'));
        load(fullfile(params.outputPath,currFieldFolder,'objects.mat'));

        processedFramesForAnalysis.params = params;
        processedFramesForAnalysis.params.inputImgsPath = fullfile(params.parentFolderForAnalysis,'Raw Images',currFieldFolder);
        processedFramesForAnalysis.objects = objects;
        processedFramesForAnalysis.frameSize = processedFrames.frameSize;
        
        %Save the objects
        save(fullfile(params.outputPath,currFieldFolder,'processedFramesForAnalysis.mat'),'processedFramesForAnalysis');  
        
        fprintf('\n%s',['SAVING OBJECT AND PROCESSED DATA FOR FIELD ',num2str(fieldNum),' of ',num2str(numel(fieldFolders)),' COMPLETE.']);

end
