function colocalMeasure = condColocGeneral(refInfo,tarInfo,condInfo,cellMask,...
                            distThreshToColocWithCond, colocDistThresh,numTarRefRand,numCondRand)
%CONDCOLOCGENERAL measure colocalization among three channels with segmented and/or detected objects 
%
%SYNOPSIS colocalMeasure = condColocGeneral(refInfo,tarInfo,condInfo,mask,...
%                          distThreshToColocWithCond, colocDistThresh,numTarRefRand,numCondRand)
%
%
%   refInfo: strcutre containing coordinates or binary image with
%            segmenation of reference channel objects
%
%   tarInfo: strcutre containing coordinates or binary image with
%            segmenation of target channel objects
%
%   condInfo: strcutre containing coordinates or binary image with
%            segmenation of condition channel objects
%
%   cellMask: binary mask of ROI (must be all ones if no real mask is input)
%             Needed to determine analyis area
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
%   numTarRefRand: number of target randomizations to assess conditional
%                  colocalization between target and reference
%                  Default 1:
%
%   numCondRand: number of condition randomizations to calculate randC
%                (see vega-Lugo et al. 2022 for randC defintion)
%                Default: 1
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
%       number of objects needed is 5)
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
%% Parse inputs

% check if reference is on the correct format
if ~isstruct(refInfo) && ~islogical(refInfo)
    error('Invalid format. refInfo must be a structure (for punctate objects) or a binary image of segmentation (for non-punctate objects)')
    
elseif isstruct(refInfo) && ~all(isfield(refInfo,{'xCoord','yCoord'}))
    error('refInfo must be a structure containing fields xCoord and yCoord')
end

if ~isstruct(tarInfo) && ~islogical(tarInfo)
    error('Invalid format. tarInfo must be a structure (for punctate objects) or a binary image of segmentation (for non-punctate objects)')
    
elseif isstruct(tarInfo) && ~all(isfield(tarInfo,{'xCoord','yCoord'}))
    error('tarInfo must be a structure containing fields xCoord and yCoord')
end

if ~isstruct(condInfo) && ~islogical(condInfo)
    error('Invalid format. condInfo must be a structure (for punctate objects) or a binary image of segmentation (for non-punctate objects)')
    
elseif isstruct(condInfo) && ~all(isfield(condInfo,{'xCoord','yCoord'}))
    error('condInfo must be a structure containing fields xCoord and yCoord')
end

%make sure a mask is input and in the correct format
if isempty(cellMask) || (~islogical(cellMask) && ~isa(cellMask,'double'))
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

if isempty(numTarRefRand)
    numTarRefRand = 1;
end

if isempty(numCondRand)
    numCondRand = 1;
end

%% Getting coordinates for each channel
imgSize = size(cellMask);

if islogical(refInfo)
    refInfo = refInfo.*cellMask;
    refSegIdx = regionprops(logical(refInfo),'PixelList','Centroid');
    refSegIdx = struct2cell(refSegIdx)';
    nRefObjects = length(refSegIdx);
    %refIsBlob = 1;
else 
    
    %convert coordinates to indices. NOTE: sub2ind assumes matrix coordinates.
    %Detection coordinates are in image coord system so we have to invert them
    %to get the right indices
    idxs = sub2ind(size(cellMask),round(refInfo.yCoord(:,1)),round(refInfo.xCoord(:,1)));

    %indices inside the mask will return 1 and 0 otherwise
    idxInside = cellMask(idxs);

    %store coords inside mask (image coords)
    refCoordVec = [refInfo.xCoord(idxInside,1) refInfo.yCoord(idxInside,1)];
    refCoordIdx = sub2ind(imgSize,round(refCoordVec(:,2)),round(refCoordVec(:,1)));
    refPtImg = zeros(imgSize);
    refPtImg(refCoordIdx) = 1;
    
    refInfo = refPtImg;
    refSegIdx = regionprops(logical(refInfo),'PixelList','Centroid');
    refSegIdx = struct2cell(refSegIdx)';
    nRefObjects = length(refSegIdx);
%     refIsBlob = 0;
end

%%%%%TAR
if islogical(tarInfo)
    tarInfo = tarInfo.*cellMask;
    tarSegIdx = regionprops(logical(tarInfo),'PixelList','Centroid');
    tarSegIdx = struct2cell(tarSegIdx)';
    nTarObjects = length(tarSegIdx);
    %tarIsBlob = 1;
else 
    
    %convert coordinates to indices. NOTE: sub2ind assumes matrix coordinates.
    %Detection coordinates are in image coord system so we have to invert them
    %to get the right indices
    idxs = sub2ind(size(cellMask),round(tarInfo.yCoord(:,1)),round(tarInfo.xCoord(:,1)));

    %indices inside the mask will return 1 and 0 otherwise
    idxInside = cellMask(idxs);

    %store coords inside mask (image coords)
    tarCoordVec = [tarInfo.xCoord(idxInside,1) tarInfo.yCoord(idxInside,1)];
    tarCoordIdx = sub2ind(imgSize,round(tarCoordVec(:,2)),round(tarCoordVec(:,1)));
    tarPtImg = zeros(imgSize);
    tarPtImg(tarCoordIdx) = 1;
    
    tarInfo = tarPtImg;
    tarSegIdx = regionprops(logical(tarInfo),'PixelList','Centroid');
    tarSegIdx = struct2cell(tarSegIdx)';
    nTarObjects = length(tarSegIdx);
    %tarIsBlob = 0;
end

if islogical(condInfo)
    condSegIdx = regionprops(logical(condInfo),'PixelList','Centroid');
    condSegIdx = struct2cell(condSegIdx)';
    nCondObjects = length(condSegIdx);
    
else 
    
    %convert coordinates to indices. NOTE: sub2ind assumes matrix coordinates.
    %Detection coordinates are in image coord system so we have to invert them
    %to get the right indices
    idxs = sub2ind(size(cellMask),round(condInfo.yCoord(:,1)),round(condInfo.xCoord(:,1)));

    %indices inside the mask will return 1 and 0 otherwise
    idxInside = cellMask(idxs);

    %store coords inside mask (image coords)
    condSegIdx = [condInfo.xCoord(idxInside,1) condInfo.yCoord(idxInside,1)];
    condCoordIdx = sub2ind(imgSize,round(condSegIdx(:,2)),round(condSegIdx(:,1)));
    condPtImg = zeros(imgSize);
    condPtImg(condCoordIdx) = 1;
    
    condInfo = condPtImg;
    
    condSegIdx = regionprops(logical(condInfo),'PixelList','Centroid');
    condSegIdx = struct2cell(condSegIdx)';
    nCondObjects = length(condSegIdx);
end

%% Create random coordinates for the condition

%when condition is blob channel
%if islogical(condInfo)
    
     realRandCondSegIdx = cell(numCondRand+1);
     realRandCondSegIdx{1} = condSegIdx;
    
    blobsizes = cellfun(@length,condSegIdx(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this was use to count blobs on condition
%flag was to tell which format to use for the condition onjects
%     numOfObjects.AllCond = nCondObjects;
%     numOfObjects.CondAreas = blobsizes;
%     flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    [~, sortIdx] = sort(blobsizes,'descend');
 for repeat = 2:numCondRand+1   
    %create the mask to be erode on each iteration(this keeps the input mask untouched)
    erodeMask = padarray(cellMask,[1 1]);
    %creat mask to store each new blob position. This mask will keep trak
    %of all new blob positions
    cellMaskWithblobs = cellMask; 
    
    %preaollocate space for random blob coordinates
    condRandCoordVec = cell(nCondObjects,1);
    blobRandShift = NaN(nCondObjects,2);

    
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
        blobRandShift(i,2) = blobRandCenter(1) - round(condSegIdx{i,1}(2));
        blobRandShift(i,1) = blobRandCenter(2) - round(condSegIdx{i,1}(1));
        
        %shift each blob pixel by blobRandShift (image coord)
        condRandCoordVec{i}(:,1) = condSegIdx{i,2}(:,1) + blobRandShift(i,1);
        condRandCoordVec{i}(:,2) = condSegIdx{i,2}(:,2) + blobRandShift(i,2);
        
        %make sure all pixels are with the mask and image limits
        condRandCoordVec{i} = max(condRandCoordVec{i},ones(size(condRandCoordVec{i}))); 
        
        condRandCoordVec{i}(:,1) = min(condRandCoordVec{i}(:,1),size(cellMask,2).*ones(size(condRandCoordVec{i}(:,1))));
        condRandCoordVec{i}(:,2) = min(condRandCoordVec{i}(:,2),size(cellMask,1).*ones(size(condRandCoordVec{i}(:,2))));
        
        %get index of shifted blob coordinates
        shiftBlobCoord = sub2ind(size(cellMask),condRandCoordVec{i}(:,2),condRandCoordVec{i}(:,1));
        %make zero the are where new blob is
        cellMaskWithblobs(shiftBlobCoord) = 0;
        erodeMask = cellMaskWithblobs;
    end
    randCondImg = ~erodeMask;
    randCondImg = randCondImg.*cellMask;
    
%     realRandCondImg(:,:,1) = condInfo;
%     realRandCondImg(:,:,2) = randCondImg;
    
    randCondSegIdx = regionprops(logical(randCondImg),'PixelList','Centroid');
    randCondSegIdx = struct2cell(randCondSegIdx)';
    
    
    realRandCondSegIdx{repeat} = randCondSegIdx;
 end    

%% Separate detections in linked and not linked

%initialize matrix for storing colocalization values. See colocalization section 
%colocalMeasure = NaN(7,2,numCondRand+1);
cmR1C1 = NaN(numCondRand+1,1);
cmR1C2 = NaN(numCondRand+1,1);
cmR2C1 = NaN(numCondRand+1,1);
cmR2C2 = NaN(numCondRand+1,1);
cmR3C1 = NaN(numCondRand+1,1);
cmR3C2 = NaN(numCondRand+1,1);
cmR4C1 = NaN(numCondRand+1,1);
cmR4C2 = NaN(numCondRand+1,1);
cmR5C1 = NaN(numCondRand+1,1);
cmR5C2 = NaN(numCondRand+1,1);
cmRow6 = NaN(numCondRand+1,1);
cmRow7 = NaN(numCondRand+1,1);

refThresh = distThreshToColocWithCond(1);
tarThresh = distThreshToColocWithCond(2);

%tempRefSegIdx = refSegIdx(:,2);

parfor i = 1:numCondRand+1
    condCoords = vertcat(realRandCondSegIdx{i}{:,2});
    
    posRefImg = zeros(imgSize);
    negRefImg = zeros(imgSize);
    
    numRefWcond = 0;
    numRefNcond = 0;
    
    for refBlob = 1:nRefObjects

       [~,d] = knnsearch(condCoords,refSegIdx{refBlob,2});    
       nnDist = mean(d);
       
       if nnDist <= refThresh
          idx = sub2ind(imgSize,refSegIdx{refBlob,2}(:,2),refSegIdx{refBlob,2}(:,1));
          posRefImg(idx) = 1;
          numRefWcond = numRefWcond + 1;
       else
           idx = sub2ind(imgSize,refSegIdx{refBlob,2}(:,2),refSegIdx{refBlob,2}(:,1));
          negRefImg(idx) = 1;
          numRefNcond = numRefNcond + 1;
       end
    end

    posTarImg = zeros(imgSize);
    negTarImg = zeros(imgSize);
    
    numTarWcond = 0;
    numTarNcond = 0;
    for tarBlob = 1:nTarObjects

       [~,d] = knnsearch(condCoords,tarSegIdx{tarBlob,2});
       nnDist = mean(d); 
       
       if nnDist <= tarThresh
          idx = sub2ind(imgSize,tarSegIdx{tarBlob,2}(:,2),tarSegIdx{tarBlob,2}(:,1));
          posTarImg(idx) = 1;
          numTarWcond = numTarWcond + 1;
       else
           idx = sub2ind(imgSize,tarSegIdx{tarBlob,2}(:,2),tarSegIdx{tarBlob,2}(:,1));
          negTarImg(idx) = 1;
          numTarNcond = numTarNcond + 1;
       end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of objects***Fix this later

     pRwC = numRefWcond/nRefObjects;
     pTwC = numTarWcond/nTarObjects;
     
     %colocalMeasure(6,1,i) = pRwC;
     %colocalMeasure(7,1,i) = pTwC;
     cmRow6(i) = pRwC;
     cmRow7(i) = pTwC;
% 
%     numOfObjects.AllTar(1,i) = totalTarDetect;
%     numOfObjects.TwC(1,i) = numTarWcondDetect;
%     numOfObjects.TnC(1,i) = numTarNcondDetect;
%     
%    numOfObjects.ROIarea = numel(find(cellMask));
%     
%     if flag
%         numOfObjects.AllRef = nBlobs;
%         numOfObjects.RefAreas = blobsizes;
%         
%     else
%         numOfObjects.AllRef(1,i) = totalRefDetect;
%         numOfObjects.RwC(1,i) = numRefWcondDetect;
%         numOfObjects.RnC(1,i) = numRefNcondDetect;
%    
%     end   
    
    %% Colocalization

    %colocalization for PosAll
    if numTarWcond >= 5 && nRefObjects >= 5
        
    [cmR1C1(i), cmR1C2(i)] = ...%[colocalMeasure(1,1,i), colocalMeasure(1,2,i)] = ...
        colocalMeasureBlob2Blob(refInfo, posTarImg,cellMask,colocDistThresh,numTarRefRand);
     
    %calculate CP/null CP
    %colocalMeasure(1,4,i) = colocalMeasure(1,1,i)/colocalMeasure(1,2,i);
    end


    %colocalization for NegAll
    if numTarNcond >= 5 && nRefObjects >= 5
        
    [cmR2C1(i), cmR2C2(i)] = ...%[colocalMeasure(2,1,i), colocalMeasure(2,2,i)] = ...
        colocalMeasureBlob2Blob(refInfo, negTarImg, cellMask,colocDistThresh,numTarRefRand);

    %calculate CP/null CP
    %colocalMeasure(2,4,i) = colocalMeasure(2,1,i)/colocalMeasure(2,2,i);
    end


    %colocalization for AllPos
    if nTarObjects >= 5 && numRefWcond >= 5
        
    [tempcmR3C1, tempcmR3C2] = ...%[colocalMeasure(3,1,i), colocalMeasure(3,2,i)] = ...
        colocalMeasureBlob2Blob(posRefImg, tarInfo, cellMask,colocDistThresh,numTarRefRand);
    
     cmR3C1(i) = tempcmR3C1/pRwC;%colocalMeasure(3,1,i) = colocalMeasure(3,1,i)/pRwC;
     cmR3C2(i) = tempcmR3C2/pRwC;%colocalMeasure(3,2,i) = colocalMeasure(3,2,i)/pRwC;
     
    %calculate CP/null CP
    %colocalMeasure(3,4,i) = colocalMeasure(3,1,i)/colocalMeasure(3,2,i);
    end


    %colocalization for AllNeg
    if nTarObjects >= 5 && numRefNcond >= 5
        
    [tempcmR4C1, tempcmR4C2] = ...%[colocalMeasure(4,1,i), colocalMeasure(4,2,i)] = ...
        colocalMeasureBlob2Blob(negRefImg, tarInfo, cellMask,colocDistThresh,numTarRefRand);
    
     cmR4C1(i) = tempcmR4C1/(1-pRwC);%colocalMeasure(4,1,i) = colocalMeasure(4,1,i)/(1-pRwC);
     cmR4C2(i) = tempcmR4C2/(1-pRwC);%colocalMeasure(4,2,i) = colocalMeasure(4,2,i)/(1-pRwC);
     
    %calculate CP/null CP
    %colocalMeasure(6,4,i) = colocalMeasure(6,1,i)/colocalMeasure(6,2,i);
    end

    % colocalization for AllAll
    if nTarObjects >= 5 && nRefObjects >= 5
        
    [cmR5C1(i), cmR5C2(i)] = ...%[colocalMeasure(5,1,i), colocalMeasure(5,2,i)] = ...
        colocalMeasureBlob2Blob(refInfo, tarInfo, cellMask,colocDistThresh,numTarRefRand);

    %calculate CP/null CP
    %colocalMeasure(7,4,i) = colocalMeasure(7,1,i)/colocalMeasure(7,2,i);
    end

end

%average coloc of all randomizations
% colocalMeasure(:,:,2) = mean(colocalMeasure(:,:,2:end),3,'omitnan');
% colocalMeasure(:,:,3:end) = [];

colocalMeasure(1,:) = [cmR1C1(1),cmR1C2(1)];
colocalMeasure(2,:) = [cmR2C1(1),cmR2C2(1)];
colocalMeasure(3,:) = [cmR3C1(1),cmR3C2(1)];
colocalMeasure(4,:) = [cmR4C1(1),cmR4C2(1)];
colocalMeasure(5,:) = [cmR5C1(1),cmR5C2(1)];
colocalMeasure(6,:) = cmRow6(1);
colocalMeasure(7,:) = cmRow7(1);

colocalMeasure(1,:,2) = mean([cmR1C1(2:end,1) cmR1C2(2:end,1)]);
colocalMeasure(2,:,2) = mean([cmR2C1(2:end,1) cmR2C2(2:end,1)]);
colocalMeasure(3,:,2) = mean([cmR3C1(2:end,1) cmR3C2(2:end,1)]);
colocalMeasure(4,:,2) = mean([cmR4C1(2:end,1) cmR4C2(2:end,1)]);
colocalMeasure(5,:,2) = mean([cmR5C1(2:end,1) cmR5C2(2:end,1)]);
colocalMeasure(6,:,2) = mean(cmRow6(2:end,1));
colocalMeasure(7,:,2) = mean(cmRow7(2:end,1));
end