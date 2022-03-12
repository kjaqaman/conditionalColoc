function condColocGeneralWrapper(MD,paramsIn)
%CONDCOLOCGENERALWRAPPER applies conditional colocalization described on condColocGeneral to a MovieData file
% movieData = condColocPt2Pt2BlobWrapper(movieData, paramsIn)
%
% Applies colocalization method described in condColocGeneral to three 
% channels. All channels can cotain either punctate or non-punctate objects
%results are stored in MovieData file directory
%
% Input:
%   movieData- A MovieData object describing the movie to be processed
%
%   paramsIn- Structure with inputs for required and optional parameters.
%   The parameters should be stored as fields in the structure, with the field
%   names and possible values as described in
%   condColocGeneralAnalysisMLMD
%       
% Output: See core function, condColocGeneral for specific outputs. 
% Outputs are saved in unique folder ConditionalColocGeneral
%   
%Jesus Vega-Lugo December 2021
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
if nargin < 1 || ~isa(MD,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end
%Get the indices of any previous colocalization processes from this function                                                                              
iProc = MD.getProcessIndex('CondColocGeneralProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(MD.processes_)+1;
    MD.addProcess(CondColocPt2Pt2BlobProcess(MD,MD.outputDirectory_));                                                                                                 
end


%Parse input, store in parameter structure
p = parseProcessParams(MD.processes_{iProc},paramsIn);

nChan = numel(MD.channels_);

if nChan ~= 3
    error('Image must have 3 channels')
end

%% Get Mask
checkMask = 1;
try
    warning('off', 'lccb:process')
    iM = MD.getProcessIndex('MaskProcess',1,0);
    channelMask = MD.getProcess(iM).funParams_.ChannelIndex;
    inMaskDir = MD.processes_{iM}.outFilePaths_(channelMask); 
    maskNames = MD.processes_{iM}.getOutMaskFileNames(channelMask);
    warning('on', 'lccb:process')
catch
    try
        warning('off', 'lccb:process')
        iM = MD.getProcessIndex('MultiThreshProcess',1,0);
        channelMask = MD.getProcess(iM).funParams_.ChannelIndex;
        inMaskDir = MD.processes_{iM}.outFilePaths_(channelMask); 
        maskNames = MD.processes_{iM}.getOutMaskFileNames(channelMask);
        warning('on', 'lccb:process')
    catch
        try
            % Try to use imported cell mask if no MaskProcess, Kevin Nguyen 7/2016
            iM = MD.getProcessIndex('ImportCellMaskProcess',Inf,0); 
            channelMask = MD.getProcess(iM).funParams_.ChannelIndex;
            inMaskDir = MD.processes_{iM}.outFilePaths_{channelMask}; 
            inMaskDir = fileparts(inMaskDir);
            inMaskDir = {inMaskDir}; % Below requires a cell
            maskNames = {{['cellMask_channel_',num2str(channelMask),'.tif']}};
        catch
            checkMask = 0;
            warning('No mask found. Using full image.')
        end
    end
end    

%% Get channel info

channelInfo = cell(1,nChan);

for ch = 1:nChan
  
    if p.PtOrBlob(ch) == 0 %is blob
        try
            iS = MD.getProcessIndex('SegmentBlobsPlusLocMaxSimpleProcess',Inf,0);
            inBlobSegDir = MD.processes_{iS}.outFilePaths_(ch); 
            load(inBlobSegDir{1},'maskBlobs');
            channelInfo{1,ch} = maskBlobs;
        catch 
            error(['No segmentation data found for channel ' num2str(ch) '! Segmentation process must be run first!'])
        end
        
    elseif p.PtOrBlob(ch) == 1%is point
        try 
           
            switch p.DetectionProcessID{ch}
               case {'SubRes' 'Subres' 'subRes' 'subres'}
                   iD = MD.getProcessIndex('SubResolutionProcess',Inf,0);
                   detectDir = MD.processes_{iD}.outFilePaths_(ch);
                   load(detectDir{1},'movieInfo');
                   channelInfo{1,ch} = movieInfo;
                   
               case {'PointSource' 'Pointsource' 'pointSource' 'pointsource'}
                   iD = MD.getProcessIndex('PointSourceDetectionProcess',Inf,0);
                   detectDir = MD.processes_{iD}.outFilePaths_(ch);
                   load(detectDir{1},'movieInfo');
                   channelInfo{1,ch} = movieInfo;
            end
           
        catch
            error(['No detection process found for channel ' num2str(ch) '! Detection process must be run first!'])
        end

    end
end

%% Set inFile and outFile Paths

condColocGeneralProc = MD.processes_{iProc};

% Set up the input directories
inFilePaths = cell(1,nChan);

for i = 1:nChan
    
    %inFilePath
    inFilePaths{1,i} = MD.getChannelPaths{i}; 
    condColocGeneralProc.setInFilePaths(inFilePaths);
end

%outFilePath
currDir = [p.OutputDirectory '/ConditionalColocGeneral'];
outFilePath = fullfile(currDir,'colocalInfo.mat');
condColocGeneralProc.setOutFilePaths(outFilePath)

mkClrDir(currDir);

%% Run Colocalization Analysis

%Load the mask for this frame/channel
if checkMask ~= 0
    currMask = imread([inMaskDir{1} filesep maskNames{1}{1}]);
else
    currMask = ones(MD.imSize);
end
currMask = logical(currMask);


%Run Function
%condColocGeneral(Ref,Tar,Cond...
colocalMeasureCond3Ref1Tar2 = condColocGeneral(channelInfo{1},...
    channelInfo{2}, channelInfo{3},currMask, p.ColocDistThresh([2 3]),p.ColocDistThresh(1),p.NumTarRefRand,p.NumCondRand);


%flip Ref and Tar. Now detectionTar is Ref and detectionRef is Tar
colocalMeasureCond3Ref2Tar1 = condColocGeneral(channelInfo{2},...
    channelInfo{1},channelInfo{3},currMask, p.ColocDistThresh([3 2]),p.ColocDistThresh(1),p.NumTarRefRand,p.NumCondRand);


colocalMeasureCond1Ref3Tar2 = condColocGeneral(channelInfo{3},...
    channelInfo{2},channelInfo{1},currMask, p.ColocDistThresh([2 1]),p.ColocDistThresh(3),p.NumTarRefRand,p.NumCondRand);


colocalMeasureCond1Ref2Tar3 = condColocGeneral(channelInfo{2},...
    channelInfo{3},channelInfo{1},currMask, p.ColocDistThresh([1 2]),p.ColocDistThresh(3),p.NumTarRefRand,p.NumCondRand);


colocalMeasureCond2Ref3Tar1 = condColocGeneral(channelInfo{3},...
    channelInfo{1},channelInfo{2},currMask, p.ColocDistThresh([3 1]),p.ColocDistThresh(2),p.NumTarRefRand,p.NumCondRand);


colocalMeasureCond2Ref1Tar3 = condColocGeneral(channelInfo{1},...
    channelInfo{3},channelInfo{2},currMask, p.ColocDistThresh([1 3]),p.ColocDistThresh(2),p.NumTarRefRand,p.NumCondRand);

save(outFilePath,'colocalMeasureCond2Ref3Tar1','colocalMeasureCond2Ref1Tar3',...
    'colocalMeasureCond3Ref2Tar1','colocalMeasureCond3Ref1Tar2',...
    'colocalMeasureCond1Ref2Tar3','colocalMeasureCond1Ref3Tar2');


disp('Finished condColocGeneral Analysis!')
end