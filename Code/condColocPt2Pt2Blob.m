function [colocalMeasure, numOfObjects] = condColocPt2Pt2Blob(refCoords,tarCoords,condCoords,...
                            mask,distThreshToColocWithCond, colocDistThresh,numRandomizations, alpha)
%CONDCOLOCPT2PT2BLOB runs conditional colocalization for 2 channels of punctate objects with one channel non-punctate object 
%
%SYNOPSIS  [colocalMeasure, numOfObjects] = condColocPt2Pt2Blob(refCoords,tarCoords,condCoords,...
%                            mask,distThreshToColocWithCond, colocDistThresh, alpha)
%
% Function divides the reference and target population into colocalized or
% not with condition. Then, calculates colocalization measures for different
% combinations of target and reference with and without condition. One
% randomization is made by randomly replacing non-punctate objects within  
% the cell area (when non-punctate object is condition) or by taking a sample  
% of points within the cell area (when a punctate object is condition).
%    
%INPUT
%   refCoords:      structure containing coordinates from reference channel
%                   If reference is non-punctate, input must be a cell
%                   array with:  
%
%                           column 1: containing coordinates of
%                           non-punctate objects
%
%                           column 2: containing centroids of non-punctate
%                           objects 
%
%                           rows: number of non-punctate objects
%
%   tarCoords:      structure containing coordinates from target channel 
%
%   condCoords:     structure containing coordinates from condition channel 
%                   If condition is non-punctate, input must be a cell
%                   array with:  
%
%                           column 1: contains coordinates of
%                           non-punctate objects
%
%                           column 2: contains centroids of non-punctate
%                           objects 
%
%                           rows: number of non-punctate objects
%
%                   Punctate object coordinates should be in image coord  
%                   system, following format of movieInfo output of
%                   detectSubResFeatures2D_StandAlone.
%
%   mask:           binary mask of ROI (must be all 1's if no real mask).
%                   Needed to determine analysis area
%
%distThreshToColocWithCond: Vector of distance threshold 
%                           for colocalization with condition.
%
%                           First entry is for reference.
%                           Second entry is for target. 
%
%                           If empty, both entries will be equal to
%                           colocDistThresh.
%
%colocDistThresh:    Distance threshold for target-reference colocalization                            
%                   (see function colocalMeasurePt2Pt)
%
%Optional
%   numRandomizations: number of randomizations to be used for calculating
%                      randC value. (see vega-Lugo et al. 2022 for randC defintion)
%                      Default: 100
%
%   alpha:          significance value for colocalization
%                   (see function colocalMeasurePt2Pt)
%                   Default: 0.05
%                   NOTE: for conditional colocalization analysis
%                         as described in Vega-Lugo et al. 2022 this
%                         parameter is not relevant.
%
%OUTPUT
%   colocalMeasure: Three dimensional matrix containing various conditional
%                   colocalization measures.
%
%                   Rows contain the following (please see Vega-Lugo et al.
%                   2022 for more detailed expalnation of below measures):
%
%                         Row1: p(TwR|TwC): probability of target
%                               colocalizing with reference given that 
%                               target colocalizes with condition.
%
%                         Row2: p(TwR|TnC): probability of target
%                               colocalizing with reference given that
%                               target does not colocalize with condition.
%                         
%                         Row3: p^rs(Tw(RwC)): rescaled probability of 
%                               target colocalizing with a reference that
%                               is itself colocalizing with condition.
%
%                         Row4: p^rs(Tw(RnC)): rescaled probability of 
%                               target colocalizing with a reference that
%                               is itself not colocalizing with condition.
%
%                         Row5: p(TwR): probability of target colocalizing
%                               with reference (regardless of either's
%                               colocalization with condition)
%
%                         Row6: p(RwC): probability of reference
%                               colocalizing with condition.
%
%                         Row7: p(TwC): probability of target
%                               colocalizing with condition.
%
%                   Columns contain the following:
%                        Column 1: colocalization measure for original
%                                  positions of reference and target objects
%                        Column 2: colocalization measure after replacing
%                                  target location with grid (nullTR).
%                        Column 3: critical value of colocalization measure
%                                  given input alpha as calculated in
%                                  Helmuth et al. BMC Bioinformatics 2010.
%                                  NOTE: this value is has not been
%                                  validated for rescaled probabilities. It 
%                                  is not used for the analysis in Vega-Lugo
%                                  et al. 2022.
%                        Column 4: ratio of column one to column two  
%           
%                   Third dimesion contains:
%                       colocalMeasure(:,:,1) contains above informaiton for 
%                                             original condition object positons
%                       colocalMeasure(:,:,2) contains above information for 
%                                             after randomizing condition
%                                             objects (randC)
%
%       NOTE: Rows 6 and 7 will only contain values in column one. For all
%       other rows, NaN indicates not enough objects for analysis (minimum
%       number of objects needed is 20)
%
%   **NOTE: When reference is non-punctate, rows 3, 4, 5, 6 and 8 from
%           above output must be ignored.
%
%   numOfObjects: structure containing below fields:
%                   AllTar: total number of target objects
%                      TwC: number of targets with condition
%                      TnC: number of target not with condition
%                   AllRef: total number of reference objects
%                      RwC: number of reference with condition
%                      RnC: number of reference not with conditon
%                  AllCond: total number of condition objects
%                  ROIArea: cell area (pixels)
%       CondAreas/RefAreas: areas of non-punctate objects (pixels) 
%
%              **Target and reference with or without condition contain 
%                a vector where the first entry shows number of molecules 
%                from the data. Second entry shows average number of
%                molecules after numRandomization of the condition.
%
%Jesus Vega-Lugo August 2019
%
%Jesus Vega-Lugo January 2021 added calculation for pRwC and pTwC. Modified
%randomization of condition when it is punctate object to erode the mask
%so randomized punctate condition positions don't fall close the edge of
%the cell mask

%% Parse inputs

% check if reference is on the correct format
if ~isstruct(refCoords) && ~iscell(refCoords)
    error('Invalid format. refCoords must be a structure (for punctate objects) or a cell (for non-punctate objects)')
    
elseif isstruct(refCoords) && ~all(isfield(refCoords,{'xCoord','yCoord'}))
    error('refCoords must be a structure containing fields xCoord and yCoord')
    
elseif iscell(refCoords) && size(refCoords,2) ~= 2
    error('refCoords must be a cell array with two columns containing non-punctate objects coordinates and centroids. See documentation for details.')
end

% check if target is on the correct format
if ~isstruct(tarCoords) || ~all(isfield(tarCoords,{'xCoord','yCoord'}))
    error('tarCoords must be a structure containing fields xCoord and yCoord')
end

% check if condition is on the correct format
if ~isstruct(condCoords) && ~iscell(condCoords)
    error('Invalid format. condCoords must be a structure (for punctate objects) or a cell (for non-punctate objects)')
    
elseif isstruct(condCoords) && ~all(isfield(condCoords,{'xCoord','yCoord'}))
    error('condCoords must be a structure containing fields xCoord and yCoord')
    
elseif iscell(condCoords) && size(condCoords,2) ~= 2
    error('condCoords must be a cell array with two columns containing non-punctate objects coordinates and centroids. See documentation for details.')
end

%make sure a mask is input and in the correct format
if isempty(mask) || (~islogical(mask) && ~isa(mask,'double'))
    error('A mask defining the are for analysis must be input')
end

%when no distThreshToColocWithCond is input 
if isempty(distThreshToColocWithCond)
    distThreshToColocWithCond(1:2) = colocDistThresh; 
end

%show error when not colocDistThresh is input
if isempty(colocDistThresh)
    error('Colocalization threshold must be input')
end

%set number of randomizations
if isempty(numRandomizations)
    numRandomizations = 100;
end

%set default alpha
if isempty(alpha)
    alpha = 0.05;
end

%% Getting coordinates for each channel

%when blobs are the reference channel

if iscell(refCoords)
    isRefBlob = 1;
    
    nBlobs = length(refCoords);
    blobsizes = cellfun(@length,refCoords(:,1));
    %get reference channel coordinates inside mask
    
    %convert cell with coordinates of each blob to a matrix. By doing this
    %information of indivual blobs is lost but is not needed, from now they
    %will be trated as a list of pixels
    refCoords = cell2mat(refCoords(:,1));
    
    %convert coordinates to indicies
    idxs = sub2ind(size(mask),round(refCoords(:,2)),round(refCoords(:,1)));
    
    %get indices inside the mask
    idxInside = mask(idxs);
    
    %store reference coordinates (image coords)
    refCoordVec = [refCoords(idxInside,1),refCoords(idxInside,2)];
    

    %get condition channel coordinates
    
    %convert coordinates to indices. See note above
    idxs = sub2ind(size(mask),round(condCoords.yCoord(:,1)),round(condCoords.xCoord(:,1)));
    
    %indices inside the mask will return 1 and 0 otherwise
    idxInside = mask(idxs);
    
    %store condition coordinates (image coords)
    condCoordVec = [condCoords.xCoord(idxInside,1) condCoords.yCoord(idxInside,1)];


else %when blob channel is the condition
    isRefBlob = 0;
    %get coordinates inside the mask for referene channel

    %convert coordinates to indices. NOTE: sub2ind assumes matrix coordinates.
    %Detection coordinates are in image coord system so we have to invert them
    %to get the right indices
    idxs = sub2ind(size(mask),round(refCoords.yCoord(:,1)),round(refCoords.xCoord(:,1)));

    %indices inside the mask will return 1 and 0 otherwise
    idxInside = mask(idxs);

    %store coords inside mask (image coords)
    refCoordVec = [refCoords.xCoord(idxInside,1) refCoords.yCoord(idxInside,1)];


    %get coordinates inside the mask for condition channel

    %store coords inside mask (image coords)
    condCoordVec = condCoords;
end

%get coordinates inside the mask for target channel

%convert coordinates to indices. See note above
idxs = sub2ind(size(mask),round(tarCoords.yCoord(:,1)),round(tarCoords.xCoord(:,1)));

%indices inside the mask will return 1 and 0 otherwise
idxInside = mask(idxs);

%store coords inside mask (image coords)
tarCoordVec = [tarCoords.xCoord(idxInside,1) tarCoords.yCoord(idxInside,1)];

%store coords in format compatible with function colocalMeasurePt2Pt.
%all reference and target coords (image coords)
refCoordstruct.xCoord = refCoordVec(:,1);
refCoordstruct.yCoord = refCoordVec(:,2);
tarCoordstruct.xCoord = tarCoordVec(:,1);
tarCoordstruct.yCoord = tarCoordVec(:,2);

%% Create random coordinates for the condition

if iscell(condCoordVec)
    vec = cell2mat(condCoordVec(:,1));
    numCondCoords = length(vec);
    condRandRealCoords = NaN(numCondCoords,2,numRandomizations+1);
    condRandRealCoords(:,:,1) = cell2mat(condCoordVec(:,1));
else
    numCondCoords = length(condCoordVec);
    condRandRealCoords = NaN(numCondCoords,2,numRandomizations+1);
    condRandRealCoords(:,:,1) = condCoordVec;
end

for repeat = 2:numRandomizations+1
%when condition is blob channel
if iscell(condCoordVec)
    %number of blobs in the condition channel
    nblobs = size(condCoordVec,1);
    blobsizes = cellfun(@length,condCoordVec(:,1));
    
    numOfObjects.AllCond = nblobs;
    numOfObjects.CondAreas = blobsizes;
       
    [~, sortIdx] = sort(blobsizes,'descend');
    
    %create the mask to be erode on each iteration(this keeps the input mask untouched)
    erodeMask = padarray(mask,[1 1]);
    %creat mask to store each new blob position. This mask will keep trak
    %of all new blob positions
    cellMaskWithblobs = mask; 
    
    %preaollocate space for random blob coordinates
    condRandCoordVec = cell(nblobs,1);
    blobRandShift = NaN(nblobs,2);

    
    for i = sortIdx' 
        %get the number of pixels to erode the mask. Divided by pi to get the radius
        blobsqrtArea = round(sqrt(blobsizes(i)/pi));
        
        erodeMask(1,:) = 0;
        erodeMask(size(erodeMask,1),:) = 0;
        erodeMask(:,1) = 0;
        erodeMask(:,size(erodeMask,2)) = 0;
        
        %erode mask
        erodeMask = imerode(erodeMask,strel('disk',blobsqrtArea));
        
        %create random centroid coords inside the mask for condition channel
        maskCoords = [];
        [maskCoords(:,1), maskCoords(:,2)] = find(erodeMask);
        %blobRandCenter in matrix coordinates
        blobRandCenter= datasample(maskCoords,1,'Replace',false);

        
        %shift every blob to have its centroid where a potar blobRandCenter is
        
        %get distance from the blob centroid to the random potar (image coord)
        blobRandShift(i,2) = blobRandCenter(1) - round(condCoordVec{i,2}(2));
        blobRandShift(i,1) = blobRandCenter(2) - round(condCoordVec{i,2}(1));
        
        %shift each blob pixel by blobRandShift (image coord)
        condRandCoordVec{i}(:,1) = condCoordVec{i,1}(:,1) + blobRandShift(i,1);
        condRandCoordVec{i}(:,2) = condCoordVec{i,1}(:,2) + blobRandShift(i,2);
        
        %make sure all pixels are with the mask and image limits
        condRandCoordVec{i} = max(condRandCoordVec{i},ones(size(condRandCoordVec{i}))); 
        
        condRandCoordVec{i}(:,1) = min(condRandCoordVec{i}(:,1),size(mask,2).*ones(size(condRandCoordVec{i}(:,1))));
        condRandCoordVec{i}(:,2) = min(condRandCoordVec{i}(:,2),size(mask,1).*ones(size(condRandCoordVec{i}(:,2))));
        
        %get index of shifted blob coordinates
        shiftBlobCoord = sub2ind(size(mask),condRandCoordVec{i}(:,2),condRandCoordVec{i}(:,1));
        %make zero the are where new blob is
        cellMaskWithblobs(shiftBlobCoord) = 0;
        erodeMask = cellMaskWithblobs;
    end
    
    %store rand coordinates for blobs
    condRandRealCoords(:,:,repeat) = cell2mat(condRandCoordVec);
    
    
else %when condtion is punctate

    %number of detections in the condition channel
    nCondDetect = size(condCoordVec,1);
    
    numOfObjects.AllCond = nCondDetect;
    
    %create random detection coords inside the mask for condition channel
    erodeMask = mask;
    erodeMask(1,:) = 0;
    erodeMask(size(erodeMask,1),:) = 0;
    erodeMask(:,1) = 0;
    erodeMask(:,size(erodeMask,2)) = 0;
    
    erodeMask = imerode(erodeMask,strel('square',3));
    
    [maskCoords(:,1), maskCoords(:,2)] = find(erodeMask);

    condRandRealCoords(:,:,repeat) = datasample(maskCoords,nCondDetect,'Replace',false);

end
end
%% Separate detections in linked and not linked
numTwC = NaN(1,numRandomizations+1);
numTnC = NaN(1,numRandomizations+1);

numRwC = NaN(1,numRandomizations+1);
numRnC = NaN(1,numRandomizations+1);

totalRefDetect = size(refCoordVec,1);
totalTarDetect = size(tarCoordVec,1);

numOfObjects.AllTar(1,1) = totalTarDetect;

%initialize matrix for storing colocalization values. See colocalization section 
colocalMeasure = NaN(7,4,numRandomizations+1);

for i = 1:numRandomizations+1
    %create a distance matrix for ref and tar with the condition
    refCondDist = distMat2(refCoordVec,condRandRealCoords(:,:,i));
    tarCondDist = distMat2(tarCoordVec,condRandRealCoords(:,:,i));

    %eliminate values above the threshold
    refCondDist(refCondDist > distThreshToColocWithCond(1)) = nan;
    tarCondDist(tarCondDist > distThreshToColocWithCond(2)) = nan;

    %get linked (Pos) and not linked (Neg) reference detections
    minDist = min(refCondDist,[],2,'omitnan');
    refWcondIdx = find(~isnan(minDist));
    refNcondIdx = find(isnan(minDist));

    %get linked (Pos) and not linked (Neg) target detections
    minDist = min(tarCondDist,[],2,'omitnan');
    tarWcondIdx = find(~isnan(minDist));
    tarNcondIdx = find(isnan(minDist));

    %store all coordinates in structure format. This is to make it compatible with
    %the colocalization code (colocalMeasurePt2Pt)

    %linked and not linked refrence coords (image coords)
    refWcondCoords.xCoord = refCoordVec(refWcondIdx,1);
    refWcondCoords.yCoord = refCoordVec(refWcondIdx,2);
    refNcondCoords.xCoord = refCoordVec(refNcondIdx,1);
    refNcondCoords.yCoord = refCoordVec(refNcondIdx,2);

    %linked and not linked target coords (image coords)
    tarWcondCoords.xCoord = tarCoordVec(tarWcondIdx,1);
    tarWcondCoords.yCoord = tarCoordVec(tarWcondIdx,2);
    tarNcondCoords.xCoord = tarCoordVec(tarNcondIdx,1);
    tarNcondCoords.yCoord = tarCoordVec(tarNcondIdx,2);
    
    numRefWcondDetect = length(refWcondIdx);
    numRefNcondDetect = length(refNcondIdx);
    
    numTarWcondDetect = length(tarWcondIdx);
    numTarNcondDetect = length(tarNcondIdx);
    
    pRwC = numRefWcondDetect/totalRefDetect;
    pTwC = numTarWcondDetect/totalTarDetect;
    
    colocalMeasure(6:7,1,i) = [pRwC; pTwC];
    
    numTwC(1,i) = numTarWcondDetect;
    numTnC(1,i) = numTarNcondDetect;

    numOfObjects.ROIarea = numel(find(mask));
    
    if isRefBlob
        numOfObjects.AllRef = nBlobs;
        numOfObjects.RefAreas = blobsizes;
        
    else
        numOfObjects.AllRef(1,1) = totalRefDetect;
        numRwC(1,i) = numRefWcondDetect;
        numRnC(1,i) = numRefNcondDetect;
   
    end   
    
    %% Colocalization

    %colocalization for PosAll
    if numTarWcondDetect >= 20 && totalRefDetect >= 20
        
    [colocalMeasure(1,1,i), colocalMeasure(1,2,i), colocalMeasure(1,3,i)] = ...
        colocalMeasurePt2Pt(refCoordstruct, tarWcondCoords, colocDistThresh(1), mask, alpha);
     
    %calculate CP/null CP
    colocalMeasure(1,4,i) = colocalMeasure(1,1,i)/colocalMeasure(1,2,i);
    end


    %colocalization for NegAll
    if numTarNcondDetect >= 20 && totalRefDetect >= 20
        
    [colocalMeasure(2,1,i), colocalMeasure(2,2,i), colocalMeasure(2,3,i)] = ...
        colocalMeasurePt2Pt(refCoordstruct, tarNcondCoords, colocDistThresh(1), mask, alpha);

    %calculate CP/null CP
    colocalMeasure(2,4,i) = colocalMeasure(2,1,i)/colocalMeasure(2,2,i);
    end


    %colocalization for AllPos
    if totalTarDetect >= 20 && numRefWcondDetect >= 20 && ~isRefBlob
        
    [colocalMeasure(3,1,i), colocalMeasure(3,2,i), colocalMeasure(3,3,i)] = ...
        colocalMeasurePt2Pt(refWcondCoords, tarCoordstruct, colocDistThresh(1), mask, alpha);
    
     colocalMeasure(3,1,i) = colocalMeasure(3,1,i)/pRwC;
     colocalMeasure(3,2,i) = colocalMeasure(3,2,i)/pRwC;
     
    %calculate CP/null CP
    colocalMeasure(3,4,i) = colocalMeasure(3,1,i)/colocalMeasure(3,2,i);
    end


    %colocalization for AllNeg
    if totalTarDetect >= 20 && numRefNcondDetect >= 20 && ~isRefBlob
        
    [colocalMeasure(4,1,i), colocalMeasure(4,2,i), colocalMeasure(4,3,i)] = ...
        colocalMeasurePt2Pt(refNcondCoords, tarCoordstruct, colocDistThresh(1), mask, alpha);
    
     colocalMeasure(4,1,i) = colocalMeasure(4,1,i)/(1-pRwC);
     colocalMeasure(4,2,i) = colocalMeasure(4,2,i)/(1-pRwC);
     
    %calculate CP/null CP
    colocalMeasure(4,4,i) = colocalMeasure(4,1,i)/colocalMeasure(4,2,i);
    end

    % colocalization for AllAll
    if totalTarDetect >= 20 && totalRefDetect >= 20
        
    [colocalMeasure(5,1,i), colocalMeasure(5,2,i), colocalMeasure(5,3,i)] = ...
        colocalMeasurePt2Pt(refCoordstruct, tarCoordstruct, colocDistThresh(1), mask, alpha);

    %calculate CP/null CP
    colocalMeasure(5,4,i) = colocalMeasure(5,1,i)/colocalMeasure(5,2,i);
    end

end

%store number of objects
numOfObjects.TwC(1,1) = numTwC(1,1);
numOfObjects.TwC(1,2) = mean(numTwC(1,2:end),'omitnan');

numOfObjects.TwC(1,1) = numTnC(1,1);
numOfObjects.TnC(1,2) = mean(numTnC(1,2:end),'omitnan');

numOfObjects.RwC(1,1) = numRwC(1,1);
numOfObjects.RwC(1,2) = mean(numRwC(1,2:end),'omitnan');

numOfObjects.RnC(1,1) = numRnC(1,1);
numOfObjects.RnC(1,2) = mean(numRnC(1,2:end),'omitnan');

%average coloc of all randomizations
colocalMeasure(:,:,2) = mean(colocalMeasure(:,:,2:end),3,'omitnan');
colocalMeasure(:,:,3:end) = [];
end