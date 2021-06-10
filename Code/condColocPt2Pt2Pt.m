function [colocalMeasure, numOfObjects] = condColocPt2Pt2Pt(refCoords,tarCoords, condCoords,...
                                mask,distThreshToColocWithCond, colocDistThresh, alpha)
%CONDCOLOCPT2PT2PT runs conditional colocalization for 3 channel images with 3 punctate objects
%
%SYNOPSIS  [colocalMeasure, numOfDetections] = condColocPt2Pt2Pt(refCoords,tarCoords, condCoords,...
%                                mask,distThreshToColocWithCond, colocDistThresh, alpha)
%
% Function divides the reference and target population into colocalized or
% not with condition. Then, calculates colocalization measures for different
% combinations of target and reference with and without condition. One
% randomization of the condition objects is made by taking a sample of  
% points within the cell area and using their coordinates as new positions
% for condition objects.
%    
%INPUT
%   refCoords:      structure containing coordinates from reference channel
%
%   tarCoords:      structure containing coordinates from target channel 
%
%   condCoords:     structure containing coordinates from condition channel
%
%                   All of the above should be in image coord system, 
%                   following format of movieInfo output of
%                   detectSubResFeatures2D_StandAlone.
%
%   mask:           binary mask of ROI (must be all 1's if no real mask).
%                   Needed to determine analysis area.
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
%   alpha:          significance value for colocalization
%                   (see function colocalMeasurePt2Pt)
%                   Default: 0.05
%                   NOTE: for conditional colocalization analysis
%                         as described in Vega-Lugo et al. 2021 this
%                         parameter is not relevant.
%OUTPUT
%   colocalMeasure: Three dimensional matrix containing various conditional
%                   colocalization measures.
%
%                   Rows contain the following (please see Vega-Lugo et al.
%                   2021 for more detailed expalnation of below measures):
%    
%                         Row1: p(TwR|TwC): probability of target
%                               colocalizing with reference given that 
%                               target colocalizes with condition.
%
%                         Row2: p(TwR|TnC): probability of target
%                               colocalizing with reference given that
%                               target does not colocalize with condition.
%
%                         Row3: p^rs(Tw(RwC)|TwC): rescaled probability of 
%                               target colocalizing with a reference that
%                               is itself colocalizing with condition given 
%                               that target is colocalizing with conditon.
%
%                         Row4: p^rs(Tw(RnC)|TnC): rescaled probability of 
%                               target colocalizing with a reference that
%                               is itself not colocalizing with condition 
%                               given that target is not colocalizing with 
%                               conditon.
%
%                         Row5: p^rs(Tw(RwC)): rescaled probability of 
%                               target colocalizing with a reference that
%                               is itself colocalizing with condition.
%
%                         Row6: p^rs(Tw(RnC)): rescaled probability of 
%                               target colocalizing with a reference that
%                               is itself not colocalizing with condition.
%
%                         Row7: p(TwR): probability of target colocalizing
%                               with reference (regardless of either's
%                               colocalization with condition)
%
%                         Row8: p(RwC): probability of reference
%                               colocalizing with condition.
%
%                         Row9: p(TwC): probability of target
%                               colocalizing with condition.
%
%                   Columns contain the following:
%                        Column 1: colocalization measure for original
%                                  positions of reference and target objects
%                        Column 2: colocalization measure after replacing
%                                  target locations with grid (nullTR).
%                        Column 3: critical value of colocalization measure
%                                  given input alpha as calculated in
%                                  Helmuth et al. BMC Bioinformatics 2010.
%                                  NOTE: this value has not been validated
%                                  for rescaled probabilities and it is
%                                  not used for the analysis in Vega-Lugo
%                                  et al. 2021.
%                        Column 4: ratio of column one to column two  
%           
%                   Third dimesion contains:
%                       colocalMeasure(:,:,1) contains above informaiton for 
%                                             original condition object positons
%                       colocalMeasure(:,:,2) contains above information for 
%                                             after randomizing condition
%                                             objects (randC)
%
%       NOTE: Rows 8 and 9 will only contain values in column one. For all
%       other rows, NaN indicates not enough objects for analysis (minimum
%       number of objects needed is 20)
%                   
%   numOfDetections: structure containing the number of objects at each 
%                    channel and category. See below fields:
%                   AllTar: total number of target objects
%                      TwC: number of targets with condition
%                      TnC: number of target not with condition
%                   AllRef: total number of reference objects
%                      RwC: number of reference with condition
%                      RnC: number of reference not with conditon
%                  AllCond: total number of condition objects
%                  ROIArea: cell area (pixels)
%                  **Target and reference with or without condition contain 
%                  two entries. The first entry shows number of objects 
%                  using original condition positions. Second entry
%                  shows number of objects after randomizing the condition.
%
%Jesus Vega-Lugo June 2019
%
%Jesus Vega-Lugo January 2021 added calculation for pRC and pTC. Modified
%randomization of condition to erode the mask so randomized punctate 
%condition positions don't fall close the edge of the cell mask
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
if ~isstruct(refCoords) || ~all(isfield(refCoords,{'xCoord','yCoord'}))
    error('refCoords must be a structure containing fields xCoord and yCoord')
end

% check if target is on the correct format
if ~isstruct(tarCoords) || ~all(isfield(tarCoords,{'xCoord','yCoord'}))
    error('tarCoords must be a structure containing fields xCoord and yCoord')
end

% check if condition is on the correct format
if ~isstruct(condCoords) || ~all(isfield(condCoords,{'xCoord','yCoord'}))
    error('condCoords must be a structure containing fields xCoord and yCoord')
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

%set default alpha
if isempty(alpha)
    alpha = 0.05;
end

%% Getting coordinates for each channel

%get coordinates inside the mask for referene channel

%convert coordinates to indices. NOTE: sub2ind assumes matrix coordinates.
%Detection coordinates are in image coord system so we have to invert them
%to get the right indices
idxs = sub2ind(size(mask),round(refCoords.yCoord(:,1)),round(refCoords.xCoord(:,1)));

%indices inside the mask will return 1 and 0 otherwise
idxInside = mask(idxs);

%store coords inside mask (matrix coords)
refCoordVec = [refCoords.yCoord(idxInside,1) refCoords.xCoord(idxInside,1)];

%get coordinates inside the mask for target channel

%convert coordinates to indices. NOTE: sub2ind assumes matrix coordinates.
%Detection coordinates are in image coord system so we have to invert them
%to get the right indices
idxs = sub2ind(size(mask),round(tarCoords.yCoord(:,1)),round(tarCoords.xCoord(:,1)));

%indices inside the mask will return 1 and 0 otherwise
idxInside = mask(idxs);

%store coords inside mask (matrix coords)
tarCoordVec = [tarCoords.yCoord(idxInside,1) tarCoords.xCoord(idxInside,1)];

%get coordinates inside the mask for condition channel

%convert coordinates to indices. NOTE: sub2ind assumes matrix coordinates.
%Detection coordinates are in image coord system so we have to invert them
%to get the right indices
idxs = sub2ind(size(mask),round(condCoords.yCoord(:,1)),round(condCoords.xCoord(:,1)));

%indices inside the mask will return 1 and 0 otherwise
idxInside = mask(idxs);

%store coords inside mask (matrix coords)
condCoordVec = [condCoords.yCoord(idxInside,1) condCoords.xCoord(idxInside,1)];

%store coords in format compatible with function colocalMeasurePt2Pt.
%all reference and observed coords (image coords)
refCoordStruct.xCoord = refCoordVec(:,2);
refCoordStruct.yCoord = refCoordVec(:,1);
tarCoordStruct.xCoord = tarCoordVec(:,2);
tarCoordStruct.yCoord = tarCoordVec(:,1);

%% Create random coordinates for the condition

%number of detections in the condition channel
nCondDetect = size(condCoordVec,1);

%stores coordinates of the data 
condRandRealCoords(:,:,1) = condCoordVec;

%create random detection coords inside the mask for condition channel
erodeMask = mask;
erodeMask(1,:) = 0;
erodeMask(size(erodeMask,1),:) = 0;
erodeMask(:,1) = 0;
erodeMask(:,size(erodeMask,2)) = 0;

erodeMask = imerode(erodeMask,strel('square',max([distThreshToColocWithCond,colocDistThresh])+3));

[maskCoords(:,1), maskCoords(:,2)] = find(erodeMask);

condRandRealCoords(:,:,2) = datasample(maskCoords,nCondDetect,'Replace',false);


%% Separate detections in linked and not linked

%initialize matrix for storing colocalization values. See colocalization section 
colocalMeasure = NaN(9,4,2);

for i = 1:2
    %create a distance matrix for ref and obs with the condition
    refCondDist = distMat2(refCoordVec,condRandRealCoords(:,:,i));
    tarCondDist = distMat2(tarCoordVec,condRandRealCoords(:,:,i));

    %eliminate values above the threshold
    refCondDist(refCondDist > distThreshToColocWithCond(1)) = nan;
    tarCondDist(tarCondDist > distThreshToColocWithCond(2)) = nan;

    %get linked (Pos) and not linked (Neg) reference detections
    minDist = min(refCondDist,[],2,'omitnan');
    refWcondIdx = find(~isnan(minDist));
    refNcondIdx = find(isnan(minDist));

    %get linked (Pos) and not linked (Neg) observed detections
    minDist = min(tarCondDist,[],2,'omitnan');
    tarWcondIdx = find(~isnan(minDist));
    tarNcondIdx = find(isnan(minDist));

    %store all coordinates in structure format. This is to make it compatible with
    %the colocalization code (colocalMeasurePt2Pt)

    %linked and not linked refrence coords (image coords)
    refWcondCoords.xCoord = refCoordVec(refWcondIdx,2);
    refWcondCoords.yCoord = refCoordVec(refWcondIdx,1);
    refNcondCoords.xCoord = refCoordVec(refNcondIdx,2);
    refNcondCoords.yCoord = refCoordVec(refNcondIdx,1);

    %linked and not linked observed coords (image coords)
    tarWcondCoords.xCoord = tarCoordVec(tarWcondIdx,2);
    tarWcondCoords.yCoord = tarCoordVec(tarWcondIdx,1);
    tarNcondCoords.xCoord = tarCoordVec(tarNcondIdx,2);
    tarNcondCoords.yCoord = tarCoordVec(tarNcondIdx,1);
    
    totalRefDetect = size(refCoordVec,1);
    numRefWcondDetect = length(refWcondIdx);
    numRefNcondDetect = length(refNcondIdx);
    
    totalTarDetect = size(tarCoordVec,1);
    numTarWcondDetect = length(tarWcondIdx);
    numTarNcondDetect = length(tarNcondIdx);
    
    pRwC = numRefWcondDetect/totalRefDetect;
    pTwC = numTarWcondDetect/totalTarDetect;
    
    colocalMeasure(8:9,1,i) = [pRwC; pTwC];
    
    %number of detections on sub populations
    numOfObjects.AllRef(1,i) = totalRefDetect;
    numOfObjects.RwC(1,i) = numRefWcondDetect;
    numOfObjects.RnC(1,i) = numRefNcondDetect;
    numOfObjects.AllTar(1,i) = totalTarDetect;
    numOfObjects.TwC(1,i) = numTarWcondDetect;
    numOfObjects.TnC(1,i) = numTarNcondDetect;
    numOfObjects.AllCond = size(condCoordVec,1);
    numOfObjects.ROIarea = numel(find(mask));
    
    %% Colocalization
    
    %colocalMeasurePt2Pt(REF,TAR,...)
    
    %colocalization for p(TwR|TwC)
    if numTarWcondDetect >= 20 && totalRefDetect >= 20
    
    [colocalMeasure(1,1,i), colocalMeasure(1,2,i), colocalMeasure(1,3,i)] = ...
        colocalMeasurePt2Pt(refCoordStruct, tarWcondCoords, colocDistThresh, mask, alpha);

    %calculate CP/null CP
    colocalMeasure(1,4,i) = colocalMeasure(1,1,i)/colocalMeasure(1,2,i);
    end


    %colocalization for p(TwR|TnC)
    if numTarNcondDetect >= 20 && totalRefDetect >= 20
        
    [colocalMeasure(2,1,i), colocalMeasure(2,2,i), colocalMeasure(2,3,i)] = ...
        colocalMeasurePt2Pt(refCoordStruct, tarNcondCoords, colocDistThresh, mask, alpha);

    %calculate cT/cNull
    colocalMeasure(2,4,i) = colocalMeasure(2,1,i)/colocalMeasure(2,2,i);
    end


    %colocalization for p^rs(Tw(RwC)|TwC)
    if numTarWcondDetect >= 20 && numRefWcondDetect >= 20

    [colocalMeasure(3,1,i), colocalMeasure(3,2,i), colocalMeasure(3,3,i)] = ...
        colocalMeasurePt2Pt(refWcondCoords, tarWcondCoords, colocDistThresh, mask, alpha);
    
    colocalMeasure(3,1,i) = colocalMeasure(3,1,i)/pRwC;
    colocalMeasure(3,2,i) = colocalMeasure(3,2,i)/pRwC;
     
    %calculate CP/null CP
    colocalMeasure(3,4,i) = colocalMeasure(3,1,i)/colocalMeasure(3,2,i);
    end


    %colocalization for p^rs(Tw(RnC)|TnC)
    if numTarNcondDetect >= 20 && numRefNcondDetect >= 20

    [colocalMeasure(4,1,i), colocalMeasure(4,2,i), colocalMeasure(4,3,i)] = ...
        colocalMeasurePt2Pt(refNcondCoords, tarNcondCoords, colocDistThresh, mask, alpha);
    
    colocalMeasure(4,1,i) = colocalMeasure(4,1,i)/(1-pRwC);
    colocalMeasure(4,2,i) = colocalMeasure(4,2,i)/(1-pRwC);

    %calculate CP/null CP
    colocalMeasure(4,4,i) = colocalMeasure(4,1,i)/colocalMeasure(4,2,i);
    end


    %colocalization for p^rs(Tw(RwC))
    if totalTarDetect >= 20 &&  numRefWcondDetect >= 20

    [colocalMeasure(5,1,i), colocalMeasure(5,2,i), colocalMeasure(5,3,i)] = ...
        colocalMeasurePt2Pt(refWcondCoords, tarCoordStruct, colocDistThresh, mask, alpha);
    
    colocalMeasure(5,1,i) = colocalMeasure(5,1,i)/pRwC;
    colocalMeasure(5,2,i) = colocalMeasure(5,2,i)/pRwC;

    %calculate CP/null CP
    colocalMeasure(5,4,i) = colocalMeasure(5,1,i)/colocalMeasure(5,2,i);
    end


    %colocalization for p^rs(Tw(RnC))
    if totalTarDetect >= 20 &&  numRefNcondDetect >= 20

    [colocalMeasure(6,1,i), colocalMeasure(6,2,i), colocalMeasure(6,3,i)] = ...
        colocalMeasurePt2Pt(refNcondCoords, tarCoordStruct, colocDistThresh, mask, alpha);
    
    colocalMeasure(6,1,i) = colocalMeasure(6,1,i)/(1-pRwC);
    colocalMeasure(6,2,i) = colocalMeasure(6,2,i)/(1-pRwC);

    %calculate CP/null CP
    colocalMeasure(6,4,i) = colocalMeasure(6,1,i)/colocalMeasure(6,2,i);
    end

    % colocalization for p(TwR)
    if totalTarDetect >= 20 && totalRefDetect >= 20
        
    [colocalMeasure(7,1,i), colocalMeasure(7,2,i), colocalMeasure(7,3,i)] = ...
        colocalMeasurePt2Pt(refCoordStruct, tarCoordStruct, colocDistThresh, mask, alpha);

    %calculate CP/null CP
    colocalMeasure(7,4,i) = colocalMeasure(7,1,i)/colocalMeasure(7,2,i);
    end

end
%substitute Infs for NaN
colocalMeasure(isinf(colocalMeasure)) = nan;
end