function movieData = condColocPt2Pt2PtWrapper(movieData, paramsIn)
%CONDCOLOCPT2PT2PtWRAPPER run conditional colocalization method described in condColocPt2Pt2Pt
%
% movieData = condColocPt2Pt2PtWrapper(movieData, paramsIn)
%
% Applies conditional colocalization described in condColocPt2Pt2Pt to 3 channel
% images. It returns the colocalization measurements described in
% colocalMeasurePt2Pt
%
%
% Input:
%   movieData- A MovieData object describing the movie to be processed
%
%   paramsIn- Structure with inputs for required and optional parameters.
%   The parameters should be stored as fields in the structure, with the field
%   names and possible values as described below.
%
%   See condColocAnalysisPt2Pt2PtMLMD documentation for parameter information   
%       
%
%
% Output: See core function, condColocPt2Pt2Pt for specific outputs. 
% Outputs are saved in unique folder ColocalizationPt2Pt2Pt
%
%
% Jesus Vega-Lugo June 2019

%% Input
% Will need to take previous process outputs from detection,segmentation and masking
% processes
%Check that input object is a valid moviedata
if nargin < 1 || ~isa(movieData,'MovieData') && ~isa(movieData,'MovieList')
    error('The first input argument must be a valid MovieData  or MovieList object!')
end

if nargin < 2
    paramsIn = [];
end
%Get the indices of any previous colocalization processes from this function
iProc = movieData.getProcessIndex('CondColocPt2Pt2PtProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(CondColocPt2Pt2PtProcess(movieData,movieData.outputDirectory_));
end


%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);


%% set outFilePath

currDir = fullfile(p.OutputDirectory,'ConditionalColocPt2Pt2Pt');
outFilePath = fullfile(currDir,'colocalInfoCond.mat');

mkdir(currDir)

%% load data

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


%Load channel 1 detection data
try
    switch p.DetectionProcessID{1} 
        
        case {'SubRes' 'Subres' 'subRes' 'subres'}
            iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
            
            channelInDetProc = NaN(size(iD));
            if channelInDetProc >1
                for iProc = 1 : length(iD)
                    channelInDetProc(iProc) = movieData.processes_{iD(iProc)}.funParams_.ChannelIndex;
                end
                iD = iD(channelInDetProc == 1);
                inDetectDir1 = movieData.processes_{iD}.outFilePaths_(1);
                
            else
                inDetectDir1 = movieData.processes_{iD}.outFilePaths_(1);
            end
           
        case {'PointSource' 'Pointsource' 'pointSource' 'pointsource'}
            iD = movieData.getProcessIndex('PointSourceDetectionProcess',Inf,0);

            inDetectDir1 = movieData.processes_{iD}.outFilePaths_(1);

    end
catch
    %Insert some error about no detection file
    error('No detection data found for Ref channel! Detection process must be run first!')
end


%Load channel 2 detection data
try
    switch p.DetectionProcessID{2} 
        
        case {'SubRes' 'Subres' 'subRes' 'subres'}
            iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
            
            channelInDetProc = NaN(size(iD));
            if channelInDetProc >1
                for iProc = 1 : length(iD)
                    channelInDetProc(iProc) = movieData.processes_{iD(iProc)}.funParams_.ChannelIndex;
                    
                    iD = iD(channelInDetProc == 2);
                    inDetectDir2 = movieData.processes_{iD}.outFilePaths_(2);
                end
            else
                inDetectDir2 = movieData.processes_{iD}.outFilePaths_(2);
            end
            
            
        
        case {'PointSource' 'Pointsource' 'pointSource' 'pointsource'}
            iD = movieData.getProcessIndex('PointSourceDetectionProcess',Inf,0);

            inDetectDir2 = movieData.processes_{iD}.outFilePaths_(2);

    end
    
catch
    %Insert some error about no detection file
    error('No detection data found for Tar channel! Detection process must be run first!')
end

%Load channel 3 detection data
try
    switch p.DetectionProcessID{3} 
        
        case {'SubRes' 'Subres' 'subRes' 'subres'}
            iD = movieData.getProcessIndex('SubResolutionProcess',Inf,0);
            
            channelInDetProc = NaN(size(iD));
            if channelInDetProc >1
                for iProc = 1 : length(iD)
                    channelInDetProc(iProc) = movieData.processes_{iD(iProc)}.funParams_.ChannelIndex;
                end
                iD = iD(channelInDetProc == 3);
                inDetectDir3 = movieData.processes_{iD}.outFilePaths_(3);
                
            else
                inDetectDir3 = movieData.processes_{iD}.outFilePaths_(3);
            end
            
        case {'PointSource' 'Pointsource' 'pointSource' 'pointsource'}
            iD = movieData.getProcessIndex('PointSourceDetectionProcess',Inf,0);

            inDetectDir3 = movieData.processes_{iD}.outFilePaths_(3);
    end
catch
    %Insert some error about no detection file
    error('No detection data found for Cond channel! Detection process must be run first!')
end

%load detection coordinates of each channel
load(inDetectDir1{1});
detectionCh1 = movieInfo;

load(inDetectDir2{1});
detectionCh2 = movieInfo;

load(inDetectDir3{1});
detectionCh3 = movieInfo;


%% Run Colocalization Analysis
 
%Load the mask for this frame/channel
if checkMask ~= 0
    currMask = imread([inMaskDir{1} filesep maskNames{1}{1}]);
else
    currMask = ones(movieData.imSize_);
end
currMask = logical(currMask);

%Run Function
%condColocPt2Pt2Pt(Ref,Tar,Cond...
[colocalMeasureCond3Ref1Tar2, numDetectCond3Ref1Tar2] = condColocPt2Pt2Pt(detectionCh1,...
    detectionCh2, detectionCh3,currMask, p.ColocDistThresh([2 3]),p.ColocDistThresh(1),p.NumRandomizations, p.AlphaValue);


%flip Ref and Tar. Now detectionTar is Ref and detectionRef is Tar
[colocalMeasureCond3Ref2Tar1, numDetectCond3Ref2Tar1] = condColocPt2Pt2Pt(detectionCh2,...
    detectionCh1, detectionCh3,currMask, p.ColocDistThresh([3 2]),p.ColocDistThresh(1),p.NumRandomizations, p.AlphaValue);


[colocalMeasureCond1Ref3Tar2, numDetectCond1Ref3Tar2] = condColocPt2Pt2Pt(detectionCh3,...
    detectionCh2,detectionCh1,currMask, p.ColocDistThresh([2 1]),p.ColocDistThresh(3),p.NumRandomizations, p.AlphaValue);


[colocalMeasureCond1Ref2Tar3, numDetectCond1Ref2Tar3] = condColocPt2Pt2Pt(detectionCh2,...
    detectionCh3,detectionCh1,currMask, p.ColocDistThresh([1 2]),p.ColocDistThresh(3),p.NumRandomizations, p.AlphaValue);


[colocalMeasureCond2Ref3Tar1, numDetectCond2Ref3Tar1] = condColocPt2Pt2Pt(detectionCh3,...
    detectionCh1,detectionCh2,currMask, p.ColocDistThresh([3 1]),p.ColocDistThresh(2),p.NumRandomizations, p.AlphaValue);


[colocalMeasureCond2Ref1Tar3, numDetectCond2Ref1Tar3] = condColocPt2Pt2Pt(detectionCh1,...
    detectionCh3,detectionCh2,currMask, p.ColocDistThresh([1 3]),p.ColocDistThresh(2),p.NumRandomizations, p.AlphaValue);

save(outFilePath,'colocalMeasureCond2Ref3Tar1','colocalMeasureCond2Ref1Tar3',...
    'colocalMeasureCond3Ref2Tar1','colocalMeasureCond3Ref1Tar2',...
    'colocalMeasureCond1Ref2Tar3','colocalMeasureCond1Ref3Tar2');

save(fullfile(currDir,'numOfObjects.mat'),'numDetectCond2Ref3Tar1',...
    'numDetectCond2Ref1Tar3','numDetectCond3Ref2Tar1','numDetectCond3Ref1Tar2',...
    'numDetectCond1Ref2Tar3','numDetectCond1Ref3Tar2');
    
%% Save Results
movieData.processes_{iProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished Colocalization Analysis!')
end