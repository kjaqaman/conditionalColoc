function [colocMeasure, numOfObjects]  = simulatePt2Pt2BlobCondColoc(numRefDetect,...
                numTarDetect,segmentation,cellMask,refColocWithCondFrac,tarColocWithCondFrac,...
                tarColocWithRefAtCondFrac,tarColocWithRefNotAtCondFrac,varargin)
%SIMULATEPT2PT2BLOBCONDCOLOC simulate object positions to be used for conditional coloalization
%
% SYNOPSIS colocMeasure, numOfObjects]  = simulatePt2Pt2BlobCondColoc(numRefDetect,...
%                numTarDetect,sementation,cellMask,refColocWithCondFrac,tarColocWithCondFrac,...
%                tarColocWithRefAtCondFrac,tarColocWithRefNotAtCondFrac,varargin)
%
% Function sets up simulated positions and runs colocalization conditional
% colocalization analysis.
%
%INPUT
% Required
%   numRefDetect: number of reference objects
%
%   numTarDetect: number of traget objects
%
%   segmentation: binary image of segmented non-punctate objects
%                NOTE: a binary image containing only points can be input.
%                      When doing this analysis will be equivalent to doing 
%                      pt2pt2pt.
%
%   cellMask: binary image of ROI
%
%   refColocWithCondFrac: scalar or vector containing the fraction of 
%                         reference objects that colocalize with condition.
%
%   tarColocWithCondFrac: scalar or vector containing the fraction of target 
%                         objects that colocalize with condition.
%
%   tarColocWithRefAtCondFrac: scalar or vector containing the fraction of 
%                              target objects colocalized with reference 
%                              objects when target colocalizes with
%                              condition.
%
%   tarColocWithRefNotAtCondFrac: scalar or vector containing the fraction 
%                                 of target objects colocalized with
%                                 reference objects when target doesn't
%                                 colocalize with condition.
% Optional enter as name-value pairs
%   colocDistRange: range of distances (in pixels) target will be
%                   positioned away from reference or segmentation objects 
%                   Default: [0 1]
%
%   colocDistThresh: Vector containing distance thresholds for conditional
%                    colocalization analysis
%                    colocDistThresh(1): threshold between 
%                                        reference and target
%                    colocDistThresh(2): threshold between reference
%                                        and segmentation
%                    colocDistThresh(3): threshold between target
%                                        and segmentation
%                    Defualt: [3 3 3] (pixels)
%OUTPUT
%   colocalMeasure: see output documentation on condColcoPt2Pt2Blob.m
%
%Jesus Vega-Lugo Dec. 2020
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

%% Parse Inputs

ip = inputParser;
ip.FunctionName = 'simulatePt2Pt2BlobCondColoc';

addRequired(ip,'numRefDetect',@isscalar)
addRequired(ip,'numTarDetect', @isscalar)
addRequired(ip,'segmentation',@islogical)
addRequired(ip,'cellMask',@(x) islogical(x) || isdouble(x))
addRequired(ip,'refColocWithCondFrac', @(x) isscalar(x) || isvector(x))
addRequired(ip,'tarColocWithCondFrac', @(x) isscalar(x) || isvector(x))
addRequired(ip,'tarColocWithRefAtCondFrac', @(x) isscalar(x) || isvector(x))
addRequired(ip,'tarColocWithRefNotAtCondFrac', @(x) isscalar(x) || isvector(x))

addParameter(ip,'colocDistRange',[0 1],@isvector)
addParameter(ip,'colocDistThresh',[3 3 3],@isvector)

parse(ip,numRefDetect,numTarDetect,segmentation,cellMask,refColocWithCondFrac,...
        tarColocWithCondFrac,tarColocWithRefAtCondFrac,...
        tarColocWithRefNotAtCondFrac,varargin{:})


numRefDetect = ip.Results.numRefDetect;
numTarDetect = ip.Results.numTarDetect;
segmentation = ip.Results.segmentation;
cellMask = ip.Results.cellMask;
refColocWithCondFrac = ip.Results.refColocWithCondFrac;
tarColocWithCondFrac = ip.Results.tarColocWithCondFrac;
tarColocWithRefAtCondFrac = ip.Results.tarColocWithRefAtCondFrac;
tarColocWithRefNotAtCondFrac = ip.Results.tarColocWithRefNotAtCondFrac;
colocDistRange = ip.Results.colocDistRange;
colocDistThresh = ip.Results.colocDistThresh;

%% Get mask coordinates

%erode mask to avoid going outside cell area
erodedMask = cellMask;
erodedMask(1,:) = 0;
erodedMask(size(erodedMask,1),:) = 0;
erodedMask(:,1) = 0;
erodedMask(:,size(erodedMask,2)) = 0;

structElem = strel('square',10);
erodedMask = imerode(erodedMask,structElem);

%get coordinates of segmentation boundaries
boundBlobs = erodedMask.*segmentation;
segBoundaries = bwboundaries(boundBlobs);
if size(segBoundaries{1},1) == 2
    for i = 1:length(segBoundaries)
        segBoundaries{i}(2,:) = [];
    end
end
segBoundaries = cell2mat(segBoundaries);

dilatedBlobs = imdilate(segmentation,structElem);
erodedMask(dilatedBlobs) = 0;

[maskWithoutBlobsCoords(:,1), maskWithoutBlobsCoords(:,2)] = find(erodedMask);
                        
%%    %%%%%%%%%%%%%%%%%%%%%%Reference%%%%%%%%%%%%%%%%%%%%%%%

    %number of colocalized ref molecules
    refDetectColocWithCond = round(numRefDetect*refColocWithCondFrac);

    %random ref positions with some noise. Total random ref
    %positions = total molecules - colocalized molecules
    refNotColocWithCondPositions = datasample(maskWithoutBlobsCoords,numRefDetect-refDetectColocWithCond,'Replace',false);
    refNotColocWithCondPositions = refNotColocWithCondPositions + rand(size(refNotColocWithCondPositions,1),2);
    
    if refDetectColocWithCond
        [randSegBoundaryCoord, condColocwithRefIdx] = datasample(segBoundaries,refDetectColocWithCond,'Replace',false);
        %refColocWithCondPositions = randSegBoundaryCoord;
        refColocWithCondPositions = randSegBoundaryCoord + (0.2*rand(refDetectColocWithCond,2)+0.2);
    else
        refColocWithCondPositions = [];
    end

    allRefPositions = vertcat(refColocWithCondPositions,refNotColocWithCondPositions);
    
    refPositionstruct.xCoord = allRefPositions(:,2);
    refPositionstruct.yCoord = allRefPositions(:,1);

%%    %%%%%%%%%%%%%%%%%%%%%%%Target%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number tar molecules colocalizing with condition +- ref
tarDetectColocWithCond = round(numTarDetect*tarColocWithCondFrac);
tarDetectColocWithRefAtCond = round(tarDetectColocWithCond*tarColocWithRefAtCondFrac);

%%%%%%%%%%%%%%%%%%%%%% tar + cond + ref%%%%%%%%%%%%%%%%%%%%
if tarDetectColocWithRefAtCond %tar colocalizing with cond and ref

    if length(refColocWithCondPositions) >= tarDetectColocWithRefAtCond 
                
        tarDetectColocWithRefAtCondPositions = refColocWithCondPositions(1:tarDetectColocWithRefAtCond,:)...
            + ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefAtCond,2) + colocDistRange(1));
    
    else
        repRefColocWithCondPositions = refColocWithCondPositions;
        while  length(repRefColocWithCondPositions) < tarDetectColocWithRefAtCond
            repRefColocWithCondPositions = repmat(repRefColocWithCondPositions,2,1);
        end

        tarDetectColocWithRefAtCondPositions = repRefColocWithCondPositions(1:tarDetectColocWithRefAtCond,:)...
        + ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefAtCond,2) + colocDistRange(1));
    end
else
    tarDetectColocWithRefAtCondPositions = [];
end

% %%%%%Check%%%%
% [~,d] = knnsearch(allRefPositions,tarDetectColocWithRefAtCondPositions,'K',1);
% [~,d1] = knnsearch(segBoundaries,tarDetectColocWithRefAtCondPositions,'K',1);
% disp(['tar with ref should ' num2str(tarDetectColocWithRefAtCond) ', is ' num2str(numel(find(d<=3))) 'and tar with cond should ' num2str(tarDetectColocWithRefAtCond) ', is ' num2str(numel(find(d1<=3)))])
% %%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% tar + cond - ref%%%%%%%%%%%%%%%%%%%%

tarDetectColocWithCondNotWithRef = tarDetectColocWithCond - tarDetectColocWithRefAtCond;

%take out cond positions where ref colocalizes
segBoundariesNoRefColoc = segBoundaries;
segBoundariesNoRefColoc(condColocwithRefIdx,:) = [];

%make sure there are enough conditions not colocalizing with reference
if length(segBoundariesNoRefColoc) >= tarDetectColocWithCondNotWithRef && ~isempty(segBoundariesNoRefColoc)
    
    condNotColocWithRef = datasample(segBoundariesNoRefColoc,...
                            tarDetectColocWithCondNotWithRef,'Replace',false);
                        
elseif ~isempty(segBoundariesNoRefColoc)
    
    repSegBoundariesNoRefColoc = segBoundariesNoRefColoc;
    while length(repSegBoundariesNoRefColoc) < tarDetectColocWithCondNotWithRef
        repSegBoundariesNoRefColoc = repmat(repSegBoundariesNoRefColoc,2,1);
    end
    condNotColocWithRef = datasample(repSegBoundariesNoRefColoc,...
                            tarDetectColocWithCondNotWithRef,'Replace',false);
    segBoundariesNoRefColoc = repSegBoundariesNoRefColoc; 
    
else
    condNotColocWithRef = [];
end
        
numCondDetectNotColocWithRef = length(condNotColocWithRef);

if tarDetectColocWithCondNotWithRef && ~isempty(condNotColocWithRef)

    if numCondDetectNotColocWithRef >= tarDetectColocWithCondNotWithRef

    %cond channel colocalized molecules with some noise. colocalization 
    %distances range from values input in colocDistRange
    tarColocWithCondNotWithRefPositions = condNotColocWithRef + ...
        ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithCondNotWithRef,2) + colocDistRange(1));
    else
        repCondDetectNoRefColoc = condNotColocWithRef;
        while length(repCondDetectNoRefColoc) < tarDetectColocWithCondNotWithRef
            repCondDetectNoRefColoc = repmat(repCondDetectNoRefColoc,2,1);
        end
        tarColocWithCondNotWithRefPositions = repCondDetectNoRefColoc + ...
        ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithCondNotWithRef,2) + colocDistRange(1));
    end
    %check if target in this category colocalizes with any reference. If
    %more than 2% of tarColocWithCondNotWithRefPositions colocalize with
    %Ref replace them for another target that doesn't colocalizes with Ref
    [~, dist] = knnsearch(allRefPositions,tarColocWithCondNotWithRefPositions,'K',1);
    numTarWithRef = sum((dist <= colocDistThresh(1))==1);
    tarColocWithCondNotWithRefPositions((dist <= colocDistThresh(1)),:) = [];
    tolerance = round(tarDetectColocWithCondNotWithRef * 0.02);
    
    if numTarWithRef > tolerance
        numTarWithRefTemp = numTarWithRef;
        
        while numTarWithRefTemp > tolerance
            for i = 1:numTarWithRefTemp

                [newTar, condNotColcWithRefIdx] = datasample(segBoundariesNoRefColoc,...
                                            1,'Replace',false);

                newTar = newTar + ((colocDistRange(2)-colocDistRange(1))*...
                        rand(1,2) + colocDistRange(1));
                
                [~, dist] = knnsearch(allRefPositions,newTar,'K',1);
                
                if dist <= colocDistThresh(3)
                    
                    segBoundariesNoRefColoc(condNotColcWithRefIdx,:) = [];
                    
                else
                    numTarWithRefTemp = numTarWithRefTemp - 1;%sum((dist <= distThreshVec(1))==1);
                    tarColocWithCondNotWithRefPositions(end+1,:) = newTar;%tarColocWithCondNotWithRefPositions(end+1:end+size(newTar,1),:) = newTar;
                end
            end
        end
    
    else
        %no need to do anything
    end
    
else
    tarColocWithCondNotWithRefPositions = [];
    warning(['No target molecules colocalize with condition but not with reference'...
        'Most likey because all condition is colocalizing with reference'])
end

% %%%%%Check%%%%
% [~,d] = knnsearch(allRefPositions,tarColocWithCondNotWithRefPositions,'K',1);
% [~,d1] = knnsearch(segBoundaries,tarColocWithCondNotWithRefPositions,'K',1);
% disp(['tar with ref should 0, is ' num2str(numel(find(d<=3))) 'and tar with cond should ' num2str(tarDetectColocWithCondNotWithRef) ', is ' num2str(numel(find(d1<=3)))])
% %%%%%%%%

tarDetectNotColocWithCond = numTarDetect-tarDetectColocWithCond;

tarDetectColocWithRefNotAtCond = round(tarDetectNotColocWithCond*tarColocWithRefNotAtCondFrac);
tarDetectNotColocWithRefOrCond = tarDetectNotColocWithCond - tarDetectColocWithRefNotAtCond;


%%%%%%%%%%%%%%%%%%%%%% tar - cond + ref%%%%%%%%%%%%%%%%%%%%
if tarDetectColocWithRefNotAtCond
    if length(refNotColocWithCondPositions) >= tarDetectColocWithRefNotAtCond

    %ref channel colocalized molecules with some noise. colocalization 
    %distances range from values input in colocDistRange
    tarColocWithRefNotWithCondPositions = refNotColocWithCondPositions(1:tarDetectColocWithRefNotAtCond,:) + ...
        ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefNotAtCond,2) + colocDistRange(1));
    else
        repRefNotColocWithCond = refNotColocWithCondPositions;
        while length(repRefNotColocWithCond) < tarDetectColocWithRefNotAtCond
            repRefNotColocWithCond = repmat(repRefNotColocWithCond,2,1);
        end
        tarColocWithRefNotWithCondPositions = repRefNotColocWithCond(1:tarDetectColocWithRefNotAtCond,:) + ...
        ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefNotAtCond,2) + colocDistRange(1));
    end
else
    tarColocWithRefNotWithCondPositions = [];
end

% %%%%%Check%%%%
% [~,d] = knnsearch(allRefPositions,tarColocWithRefNotWithCondPositions,'K',1);
% [~,d1] = knnsearch(segBoundaries,tarColocWithRefNotWithCondPositions,'K',1);
% disp(['tar with ref should ' num2str(tarDetectColocWithRefNotAtCond) ', is ' num2str(numel(find(d<=3))) 'and tar with cond should 0, is ' num2str(numel(find(d1<=3)))])
% %%%%%%%%

%%%%%%%%%%%%%%%%tar - cond - ref %%%%%%%%%%%%%%%%%%%

%erode ref positions to make sure there is no random coloc
refDetectMask = zeros(size(segmentation));
refDetectIdx = sub2ind(size(segmentation),round(allRefPositions(:,1)),round(allRefPositions(:,2)));

refDetectMask(refDetectIdx) = 1;
refDetectMask = imdilate(refDetectMask,strel('square',10));

erodedMask(logical(refDetectMask)) = 0;
[maskWithoutBlobsAndRefCoords(:,1),maskWithoutBlobsAndRefCoords(:,2)] = find(erodedMask);

%random tar positions with some noise. Total random tar
%positions = total molecules - colocalized molecules
tarDetectNotColocWithRefOrCondPositions = datasample(maskWithoutBlobsAndRefCoords,...
                            tarDetectNotColocWithRefOrCond,'Replace',false);
                        
 tarDetectNotColocWithRefOrCondPositions = tarDetectNotColocWithRefOrCondPositions...
                 + rand(size(tarDetectNotColocWithRefOrCondPositions,1),2);
% %%%%%%Check%%%
% [~,d] = knnsearch(allRefPositions,tarDetectNotColocWithRefOrCondPositions,'K',1);
% [~,d1] = knnsearch(segBoundaries,tarDetectNotColocWithRefOrCondPositions,'K',1);
% disp(['tar with ref should 0, is ' num2str(numel(find(d<=3))) 'and tar with cond should 0, is ' num2str(numel(find(d1<=3)))])
% %%%%%%%%  

allTarPositions = vertcat(tarDetectColocWithRefAtCondPositions,tarColocWithCondNotWithRefPositions,...
    tarColocWithRefNotWithCondPositions,tarDetectNotColocWithRefOrCondPositions);
    
tarPositionStruct.xCoord = allTarPositions(:,2);
tarPositionStruct.yCoord = allTarPositions(:,1);

%% Run conditional colocalization

%but first get segmentation coordinates
SegIdx = regionprops(logical(segmentation),'PixelList','Centroid');
nObjects = length(SegIdx);
segmentationData = cell(nObjects,1);
for i = 1:nObjects
    segmentationData{i,1} = SegIdx(i).PixelList;
    segmentationData{i,2} = SegIdx(i).Centroid;
end    

[colocMeasure,numOfObjects] = condColocPt2Pt2Blob(refPositionstruct,...
    tarPositionStruct, segmentationData,cellMask,...
    [colocDistThresh(2) colocDistThresh(3)],colocDistThresh(1),0.05);
    
end
