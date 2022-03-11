function segmentLocMaxWrapper(movieData, varargin)
%DETECTMOVIETHRESHLOCMAX compiles detection data from movie frames

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
 if ~iscell(thresholdMethod)
     t = thresholdMethod;
     thresholdMethod = cell(1);
     thresholdMethod{1} = t;
 end
 methodValue= p.detectionParam.methodValue;
 filterNoise=p.detectionParam.filterNoise;
 filterBackground=p.detectionParam.filterBackground; 
 minSize=p.detectionParam.minSize;
 locMax = p.detectionParam.locMax;
 plotRes = p.plotRes;
%  alphaLocMax=p.detectionParam.alphaLocMax;%Not used for now...
 
 %Was masking done?
 %Define which process was masking process
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

%% --------------- Initialization ---------------%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Jesus Vega-Lugo (12/2021) modified this and below sections for them to have
%a loop to segment multiple channels in one run instead of segmenting
%channels individually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nChan =length(p.ChannelIndex);

% Set up the input directories
inFilePaths = cell(1,nChan);

% Set up the output directories
outFilePaths = cell(1,nChan);
saveResults(nChan,1)=struct();
dName = 'blob_segmentation_for_channel_';

for i = 1:nChan
    currChan = p.ChannelIndex(i);
    
    inFilePaths{1,currChan} = movieData.getChannelPaths{currChan};
    
    segmentLocMaxProc.setInFilePaths(inFilePaths);
   
    currDir = [p.OutputDirectory filesep dName num2str(currChan)];
    saveResults(currChan).dir = currDir ;
    saveResults(currChan).filename = ['Channel_' num2str(currChan) '_segmentation_result.mat'];
    %Create string for current directory
    outFilePaths{1,currChan} = [saveResults(currChan).dir filesep saveResults(currChan).filename ];
    segmentLocMaxProc.setOutFilePaths(outFilePaths{1,currChan},currChan);
    mkClrDir(currDir);
end


%% --------------- Segmentation Local maxima detection ---------------%%% 

%load cell mask
if checkMask ~= 0
    %Load the mask for this frame/channel
    mask = imread([inMaskDir{1} filesep maskNames{1}{1}]);
else
    mask= ones(movieData.imSize_);
end

for i = 1:nChan
    ch = p.ChannelIndex(i);
    disp(['Please wait, segmenting objects for channel ' num2str(ch)])
    disp(inFilePaths{1,ch});
    disp('Results will be saved under:')
    disp(outFilePaths{1,ch});
    

        %Load Image
        I = movieData.channels_(ch).loadImage(1);
        
        %Run Detection    
    [maskBlobs,labels,maskBlobsVis] = segmentBlobs_locmax(I,...
    thresholdMethod{ch},methodValue(ch),filterNoise(ch),filterBackground(ch),minSize(ch),locMax(ch),plotRes(ch),mask); 
        
    save(strcat(saveResults(ch).dir,'/',saveResults(ch).filename),'maskBlobs','labels','maskBlobsVis');
end

end