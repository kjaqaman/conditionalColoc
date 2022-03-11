function movieData = condColocPt2Pt2BlobWrapper(movieData, paramsIn)
%CONDCOLOCPT2PT2BLOBWRAPPER applies colocalization method described in condColocPt2Pt2Blob to a MovieData file
% 
% movieData = condColocPt2Pt2BlobWrapper(movieData, paramsIn)
%
% Applies colocalization method described in condColocPt2Pt2Blob to three 
% channels (punctate (2 channels) and blob(segmentable objects; 1 channel)) 
% results are stored in MovieData file directory
%
%
% Input:
%   movieData- A MovieData object describing the movie to be processed
%
%   paramsIn- Structure with inputs for required and optional parameters.
%   The parameters should be stored as fields in the structure, with the field
%   names and possible values as described in
%   condColocAnalysisPt2Pt2BlobMLMD
%       
% Output: See core function, condColocPt2Pt2Blob for specific outputs. 
% Outputs are saved in unique folder ColocalizationPt2Pt2Blob
%   
%Jesus Vega-Lugo 08/2019
%
% Copyright (C) 2021, Jaqaman Lab - UTSouthwestern 
%
% This file is part of conditionalColoc.
% 
% conditionalColoc is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% conditionalColoc is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with conditionalColoc.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
%% Input
% Will need to take previous process outputs from detection,segmentation and masking 
% processes
%Check that input object is a valid moviedata 
    if nargin < 1 || ~isa(movieData,'MovieData')
        error('The first input argument must be a valid MovieData object!')
    end

    if nargin < 2
        paramsIn = [];
    end
%Get the indices of any previous colocalization processes from this function                                                                              
    iProc = movieData.getProcessIndex('CondColocPt2Pt2BlobProcess',1,0);

%If the process doesn't exist, create it
    if isempty(iProc)
        iProc = numel(movieData.processes_)+1;
        movieData.addProcess(CondColocPt2Pt2BlobProcess(movieData,movieData.outputDirectory_));                                                                                                 
    end


%Parse input, store in parameter structure
    p = parseProcessParams(movieData.processes_{iProc},paramsIn);

    nChan = numel(movieData.channels_);

    if max(p.DetectedChannels) > nChan || min(p.DetectedChannels)<1 || ~isequal(round(p.DetectedChannels),p.DetectedChannels)
        error('Invalid channel numbers specified! Check DetectedChannel input!!')
    end

    if max(p.ChannelBlob) > nChan || min(p.ChannelBlob)<1 || ~isequal(round(p.ChannelBlob),p.ChannelBlob)
        error('Invalid channel numbers specified! Check ChannelBlob input!!')
    end
    

    %Define which process was masking process
    
        checkMask = 1;
        try
            warning('off', 'lccb:process')
            iM = movieData.getProcessIndex('MaskProcess',1,0);
            channelMask = movieData.getProcess(iM).funParams_.ChannelIndex;
            inMaskDir = movieData.processes_{iM}.outFilePaths_(channelMask); 
            maskNames = movieData.processes_{iM}.getOutMaskFileNames(channelMask);
            warning('on', 'lccb:process')
        catch
            try
                warning('off', 'lccb:process')
                iM = movieData.getProcessIndex('MultiThreshProcess',1,0);
                channelMask = movieData.getProcess(iM).funParams_.ChannelIndex;
                inMaskDir = movieData.processes_{iM}.outFilePaths_(channelMask); 
                maskNames = movieData.processes_{iM}.getOutMaskFileNames(channelMask);
                warning('on', 'lccb:process')
            catch
                try
                    % Try to use imported cell mask if no MaskProcess, Kevin Nguyen 7/2016
                    iM = movieData.getProcessIndex('ImportCellMaskProcess',Inf,0); 
                    channelMask = movieData.getProcess(iM).funParams_.ChannelIndex;
                    inMaskDir = movieData.processes_{iM}.outFilePaths_{channelMask}; 
                    inMaskDir = fileparts(inMaskDir);
                    inMaskDir = {inMaskDir}; % Below requires a cell
                    maskNames = {{['cellMask_channel_',num2str(channelMask),'.tif']}};
                catch
                    checkMask = 0;
                    warning('No mask found. Using full image.')
                end
            end
        end

    %load point channels
        try
            switch p.DetectionProcessID{1} 
        
            case {'SubRes' 'Subres' 'subRes' 'subres'}
                iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
                
                
                if numel(iD) > 1
                    
                    channelInDetProc = NaN(size(iD));
                    for iProc = 1 : length(iD)
                        channelInDetProc(iProc) = movieData.processes_{iD(iProc)}.funParams_.ChannelIndex;
                    end
                    
                    iD = iD(channelInDetProc == p.DetectedChannels(1));
                    detect1Dir = movieData.processes_{iD}.outFilePaths_(p.DetectedChannels(1));
                    
                elseif numel(iD) == 1
                   
                   detect1Dir = movieData.processes_{iD}.outFilePaths_(1);
                end

                

            case {'PointSource' 'Pointsource' 'pointSource' 'pointsource'}
                iD = movieData.getProcessIndex('PointSourceDetectionProcess',Inf,0);

                if numel(iD) > 1
                    
                    channelInDetProc = NaN(size(iD));
                    for iProc = 1 : length(iD)
                        channelInDetProc(iProc) = movieData.processes_{iD(iProc)}.funParams_.ChannelIndex;
                    end
                    
                    iD = iD(channelInDetProc == p.DetectedChannels(1));
                    detect1Dir = movieData.processes_{iD}.outFilePaths_(p.DetectedChannels(1));
                    
                elseif numel(iD) == 1
                   channelInDetProc = movieData.processes_{iD}.funParams_.ChannelIndex;
                   chIdx = channelInDetProc(channelInDetProc == p.DetectedChannels(1));
                   detect1Dir = movieData.processes_{iD}.outFilePaths_(chIdx);
                end
            end
        
        catch 
            %Insert some error about no detection file
            error(['No detection data found for channel' num2str(p.DetectedChannels(1)) '! Detection process must be run first!'])
        end
        
        %load the other point channel
        try
            switch p.DetectionProcessID{2} 
        
            case {'SubRes' 'Subres' 'subRes' 'subres'}
                iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);

                if numel(iD) > 1
                    
                    channelInDetProc = NaN(size(iD));
                    for iProc = 1 : length(iD)
                        channelInDetProc(iProc) = movieData.processes_{iD(iProc)}.funParams_.ChannelIndex;
                    end
                    
                    iD = iD(channelInDetProc == p.DetectedChannels(2));
                    detect2Dir = movieData.processes_{iD}.outFilePaths_(p.DetectedChannels(2));
                    
                elseif numel(iD) == 1
                   channelInDetProc = movieData.processes_{iD}.funParams_.ChannelIndex;
                   chIdx = channelInDetProc(channelInDetProc == p.DetectedChannels(2));
                   detect2Dir = movieData.processes_{iD}.outFilePaths_(chIdx);
                end

                

            case {'PointSource' 'Pointsource' 'pointSource' 'pointsource'}
                iD = movieData.getProcessIndex('PointSourceDetectionProcess',Inf,0);

                if numel(iD) > 1
                    
                    channelInDetProc = NaN(size(iD));
                    for iProc = 1 : length(iD)
                        channelInDetProc(iProc) = movieData.processes_{iD(iProc)}.funParams_.ChannelIndex;
                    end
                    
                    iD = iD(channelInDetProc == p.DetectedChannels(2));
                    detect2Dir = movieData.processes_{iD}.outFilePaths_(p.DetectedChannels(2));
                    
                elseif numel(iD) == 1
                   channelInDetProc = movieData.processes_{iD}.funParams_.ChannelIndex;
                   chIdx = channelInDetProc(channelInDetProc == p.DetectedChannels(2));
                   detect2Dir = movieData.processes_{iD}.outFilePaths_(chIdx);
                end
            end
        catch 
            %Insert some error about no detection file
            error(['No detection data found for channel' num2str(p.DetectedChannels(2)) '! Detection process must be run first!'])
        end
        
    %Load blob channel segmentation data
        try
            iS = movieData.getProcessIndex('SegmentBlobsPlusLocMaxSimpleProcess',Inf,0);
            inBlobSegDir = movieData.processes_{iS}.outFilePaths_(p.ChannelBlob); 
        catch 
            %Insert some error about no segmentation file
            error('No segmentation data found! Detection process must be run first!')
        end
        load(detect1Dir{1});
        detection1 = movieInfo;
        load(detect2Dir{1});
        detection2 = movieInfo;
        load(inBlobSegDir{1});

%% Run Colocalization Analysis

        %Load the mask for this frame/channel
        if checkMask ~= 0
            currMask = imread([inMaskDir{1} filesep maskNames{1}{1}]);
        else
            currMask = ones(movieData.imSize);
        end
        currMask = logical(currMask);
        
        %Get segmenation coordinates from segmentation data
        SegIdx = regionprops(logical(maskBlobs),'PixelList','Centroid');
        nObjects = length(SegIdx);
        segmentationData = cell(nObjects,1);
        for i = 1:nObjects
            segmentationData{i,1} = SegIdx(i).PixelList;
            segmentationData{i,2} = SegIdx(i).Centroid;
        end
        
        %Run Function
        [colocalMeasureCondBlobRef1Tar2, numObjCondBlobRef1Tar2] = condColocPt2Pt2Blob(detection1, detection2,...
            segmentationData,currMask,p.ColocDistThresh([2 3]),p.ColocDistThresh(1),p.NumRandomizations,p.AlphaValue);          

        [colocalMeasureCondBlobRef2Tar1, numObjCondBlobRef2Tar1] = condColocPt2Pt2Blob(detection2, detection1,...
            segmentationData,currMask,p.ColocDistThresh([3 2]),p.ColocDistThresh(1),p.NumRandomizations,p.AlphaValue);          
        
        [colocalMeasureCond2RefBlobTar1, numObjCond2RefBlobTar1] = condColocPt2Pt2Blob(segmentationData, detection1,...
            detection2,currMask,p.ColocDistThresh([3 1]),p.ColocDistThresh(2),p.NumRandomizations,p.AlphaValue);          
        
        %colocalMeasureCond2RefBlobTar1((3:4),:) = nan; 
        
        [colocalMeasureCond1RefBlobTar2, numObjCond1RefBlobTar2] = condColocPt2Pt2Blob(segmentationData, detection2,...
            detection1,currMask,p.ColocDistThresh([2 1]),p.ColocDistThresh(3),p.NumRandomizations,p.AlphaValue); 
        
        %colocalMeasureCond1RefBlobTar2((3:4),:) = nan;
        
        colocalizationParameters.detectedChannels = p.DetectedChannels;
        colocalizationParameters.blobChannel = p.ChannelBlob;
        colocalizationParameters.colocDistThresh = p.ColocDistThresh;
        
    mkdir([p.OutputDirectory '/ColocalizationPt2Pt2Blob' ])
    save([p.OutputDirectory '/ColocalizationPt2Pt2Blob/colocalInfo.mat'],...
        'colocalMeasureCondBlobRef1Tar2','colocalMeasureCondBlobRef2Tar1',...
        'colocalMeasureCond2RefBlobTar1','colocalMeasureCond1RefBlobTar2',...
        'colocalizationParameters');
    
    save([p.OutputDirectory '/ColocalizationPt2Pt2Blob/numOfObjects.mat'],...
        'numObjCondBlobRef1Tar2','numObjCondBlobRef2Tar1','numObjCond2RefBlobTar1',...
        'numObjCond1RefBlobTar2')


%% Save Results
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished Colocalization Analysis!')