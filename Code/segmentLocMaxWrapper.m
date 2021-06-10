function segmentLocMaxWrapper(movieData, varargin)
%DETECTMOVIETHRESHLOCMAX compiles detection data from movie frames
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

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('SegmentBlobsPlusLocMaxSimpleProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SegmentLocMaxProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
segmentLocMaxProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(segmentLocMaxProc,paramsIn);

%% Detection Parameters
 thresholdMethod=p.detectionParam.thresholdMethod;
 methodValue= p.detectionParam.methodValue;
 filterNoise=p.detectionParam.filterNoise;
 filterBackground=p.detectionParam.filterBackground; 
 minSize=p.detectionParam.minSize;
 locMax = p.detectionParam.locMax;
 plotRes = p.plotRes;
%  alphaLocMax=p.detectionParam.alphaLocMax;%Not used for now...
 
 %Was masking done?
 %Define which process was masking process
 if ~isempty(p.ChannelMask)
     checkMask = 1;
    try
        warning('off', 'lccb:process')
        iM = movieData.getProcessIndex('MaskProcess',1,0); 
        channelMask = movieData.getProcess(iM).funParams_.ChannelIndex;
        inMaskDir = movieData.processes_{iM}.outFilePaths_(channelMask); 
        maskNames = movieData.processes_{iM}.getOutMaskFileNames(1);
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
                warning('No mask found. Using whole image.')
                checkMask = 0;
            end
        end
    end
 else
     checkMask = 0;
 end

%% --------------- Initialization ---------------%%

nChan=numel(movieData.channels_);
% Set up the input directories
inFilePaths = cell(1,nChan);
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
segmentLocMaxProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,nChan);
saveResults(nChan,1)=struct();
dName = 'blob_segmentation_for_channel_';
for i = p.ChannelIndex   
    currDir = [p.OutputDirectory filesep dName num2str(i)];
    saveResults(i).dir = currDir ;
    saveResults(i).filename = ['Channel_' num2str(i) '_segmentation_result.mat'];
    %Create string for current directory
    outFilePaths{1,i} = [saveResults(i).dir filesep saveResults(i).filename ];
    segmentLocMaxProc.setOutFilePaths(outFilePaths{1,i},i);
    mkClrDir(currDir);
end


%% --------------- Segmentation Local maxima detection ---------------%%% 
disp('Starting segmentation and local maxima detection...')

    disp(['Please wait, segmenting objects for channel ' num2str(p.ChannelIndex)])
    disp(inFilePaths{1,p.ChannelIndex});
    disp('Results will be saved under:')
    disp(outFilePaths{1,p.ChannelIndex});
    

        %Load Image
        I = movieData.channels_(p.ChannelIndex).loadImage(1);
        if checkMask ~= 0
            %Load the mask for this frame/channel
            mask = imread([inMaskDir{1} filesep maskNames{1}{1}]);
        else
            mask= ones(movieData.imSize_);
        end
        %Run Detection    
    [maskBlobs,labels,maskBlobsVis] = segmentBlobs_locmax(I,...
    thresholdMethod,methodValue,filterNoise,filterBackground,minSize,locMax,plotRes,mask); 
        
    save(strcat(saveResults(p.ChannelIndex).dir,'/',saveResults(p.ChannelIndex).filename),'maskBlobs','labels','maskBlobsVis');


end