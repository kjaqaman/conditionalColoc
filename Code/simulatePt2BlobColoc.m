function [pTR, nullTR, cCritical] = simulatePt2BlobColoc(numTarDetect,segmentation,...
                        cellMask,colocFraction,varargin)
% SIMULATEPT2BLOBCOLOC simulate object positions to be used for colocalization
%
% SYNOPSIS 
%INPUT
%   numTarDetect: number of target objects
%
%   segmentation: binary image of segmented non-punctate objects
%                NOTE: a binary image containing only points can be input.
%                      When doing this analysis will be equivalent to doing 
%                      pt2pt.
%
%   cellMask: binary image of ROI
%
%   colocFraction: fraction of target objects colocalizing with non-punctate objects.
%
% Optional enter as name-value pairs
%   colocDistRange: range of distances (in pixels) target will be
%                   positioned away from reference or segmentation objects 
%                   Default: [0 1]
%
%   distThresh: distance threshold for colocalization analysis
%               Default: 3
%
%OUTPUT
%   See colocalMeasurePt2Pt.m documentation for output information
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

%% Parse input

ip = inputParser;
ip.FunctionName = 'simulatePt2BlobColoc';


addRequired(ip,'numTarDetect', @(x) isscalar(x) || isvector(x))
addRequired(ip,'segmentation',@islogical)
addRequired(ip,'cellMask',@(x) islogical(x) || isdouble(x))
addRequired(ip,'colocFraction',@(x) isscalar(x) || isvector(x))

addParameter(ip,'colocDistRange',[0 1],@isvector)
addParameter(ip,'distThresh',3,@isscalar)

parse(ip,numTarDetect,segmentation,cellMask,colocFraction,varargin{:})

numTarDetect = ip.Results.numTarDetect;
segmentation = ip.Results.segmentation;
cellMask = ip.Results.cellMask;
colocFraction = ip.Results.colocFraction;
colocDistRange = ip.Results.colocDistRange;
distThresh = ip.Results.distThresh;

%% Get mask coords

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

%% Target

tarColocDetect = round(numTarDetect*colocFraction);

tarColocWithBlobPositions = datasample(segBoundaries,tarColocDetect,'Replace',true);

tarColocDetectPositions = tarColocWithBlobPositions+ ((colocDistRange(2)-colocDistRange(1))*...
                            rand(tarColocDetect,2) + colocDistRange(1));
                        
tarDetectPositions = datasample(maskWithoutBlobsCoords,numTarDetect-tarColocDetect,'Replace',false);
tarDetectPositions = tarDetectPositions + rand(size(tarDetectPositions,1),2);

allTarPositions = vertcat(tarColocDetectPositions,tarDetectPositions);

tarPositionStruct.xCoord = allTarPositions(:,2);
tarPositionStruct.yCoord = allTarPositions(:,1);

segmentationCoords = regionprops(logical(segmentation),'PixelList');
segmentationCoords = struct2cell(segmentationCoords)';
segmentationCoords = cell2mat(segmentationCoords);
segmentationCoords = segmentationCoords(:,[2 1]);

%% Run pt2blob colocalization

[pTR, nullTR, cCritical] = colocalMeasurePt2Blob(segmentationCoords,...
                            tarPositionStruct,distThresh,cellMask,0.05);
end