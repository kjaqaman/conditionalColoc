function [pTR,nullTR] = simulateBlob2BlobColoc(segmentation,numTarObj,cellMask,colocFraction,varargin)
%SIMULATEBLOB2BLOBCOLOC simulate object positions to be used for blob2blob coloalization
%
% Function sets up simulated positions and runs blob2blob colocalization analysis.
%
%INPUT
% Required
%   segmentation: binary image of segmented non-punctate objects
%                NOTE: a binary image containing only points can be input.
%                      When doing this analysis will be equivalent to doing 
%                      pt2pt2pt.
%
%   numTarObj: number of traget objects.
%
%   cellMask: binary image of ROI.
%
%   colocFraction: scalar o vector containing the desired fraction(s) of
%                  target colocalization with reference.
%Optional
%   numRandomizations: number of randomizations to be used for calculating
%                      randC value. (see vega-Lugo et al. 2022 for randC defintion)
%                      Default: 1
%
%   tarBlobOrPt: 1 to create non-punctate traget object, 0 for creating
%                punctate target object. Default: 1
%
%   lineLength: length of line to be used to dilate the center of
%               non-punctate target and reference objects. Default: 7
%
%   tarDiskSize: radius for disk to dilate initial line dilation of target
%                objects. Default: 5
%
%OUTPUT
%   pTR: probability of target colocalizing with reference.
%
%   nullTR: probability of coincidental target-reference colocalization
%           (after target randomization).
%
%Jesus vega-Lugo December 2021
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
ip.FunctionName = 'simulateBlob2BlobColoc';

addRequired(ip,'segmentation',@islogical)
addRequired(ip,'numTarObj',@isscalar)
addRequired(ip,'cellMask',@(x) islogical(x) || isdouble(x))
addRequired(ip,'colocFraction', @(x) isscalar(x) || isvector(x))

addParameter(ip,'colocDistRange',[0 1],@isvector)
addParameter(ip,'distThresh',3,@isscalar)
addParameter(ip,'numRandomizations',1,@isscalar)
addParameter(ip,'tarBlobOrPt',1,@(x) x == 0 || x == 1)
addParameter(ip,'tarDiskSize',5,@isscalar)
addParameter(ip,'lineLength',7,@isscalar)
parse(ip,segmentation,numTarObj,cellMask,colocFraction,varargin{:})


numTarObj = ip.Results.numTarObj;
segmentation = ip.Results.segmentation;
cellMask = ip.Results.cellMask;
colocFraction = ip.Results.colocFraction;
colocDistRange = ip.Results.colocDistRange;
distThresh = ip.Results.distThresh;
numRandomizations = ip.Results.numRandomizations;
tarBlobOrPt = ip.Results.tarBlobOrPt;
tarDiskSize = ip.Results.tarDiskSize;
lineLength = ip.Results.lineLength;

%% Get mask coordinates

%erode mask to avoid going outside cell area
erodedMask = cellMask;
erodedMask(1,:) = 0;
erodedMask(size(erodedMask,1),:) = 0;
erodedMask(:,1) = 0;
erodedMask(:,size(erodedMask,2)) = 0;

structElem = strel('square',15);
erodedMask = imerode(erodedMask,structElem);

erodedCellMask = erodedMask;

%get coordinates of segmentation boundaries
boundBlobs = erodedMask.*segmentation;
segBoundaries = bwboundaries(boundBlobs);
if size(segBoundaries{1},1) == 2
    for i = 1:length(segBoundaries)
        segBoundaries{i}(2,:) = [];
    end
end

coordsPerSeg = cellfun(@length,segBoundaries);
segBoundaries = cell2mat(segBoundaries);

dilatedBlobs = imdilate(segmentation,structElem);
erodedMask(dilatedBlobs) = 0;

[maskWithoutBlobsCoords(:,1), maskWithoutBlobsCoords(:,2)] = find(erodedMask);
                        
%%    %%%%%%%%%%%%%%%%%%%%%%Target%%%%%%%%%%%%%%%%%%%%%%%
imgSize = size(cellMask);

tarWrefImg = zeros(imgSize);
tarNotwRefImg = zeros(imgSize);

%number of colocalized tar molecules
 tarObjColocWithRef = round(numTarObj*colocFraction);

 numTarNotWref = numTarObj-tarObjColocWithRef;
 
 tarNotColocWithRef = NaN(numTarNotWref,2);
 
 tarPtMask = erodedCellMask;
 
 for i = 1:numTarNotWref
    tarNotColocWithRef(i,:) = datasample(maskWithoutBlobsCoords,1,'Replace',false);
    maskWithoutBlobsCoords = [];
    tarPtMask(tarNotColocWithRef(i,1),tarNotColocWithRef(i,2)) = 0;
    tarMask = imerode(tarPtMask,strel('disk',tarDiskSize+5));
    tarMaskWithoutBlobs = tarMask.*erodedMask;
    [maskWithoutBlobsCoords(:,1), maskWithoutBlobsCoords(:,2)] = find(tarMaskWithoutBlobs);
 end   

 tarNotColocWithRefIdx = sub2ind(size(cellMask),tarNotColocWithRef(:,1),tarNotColocWithRef(:,2));
 
 tarNotwRefImg(tarNotColocWithRefIdx) = 1;


 if tarObjColocWithRef
     tarColocWithRef = NaN(tarObjColocWithRef,2);
     numSeg = length(coordsPerSeg); 
     
     if numSeg < tarObjColocWithRef
         coordsPerSeg = vertcat(coordsPerSeg,coordsPerSeg);
         segBoundaries = vertcat(segBoundaries,segBoundaries);
         
         numSeg = length(coordsPerSeg);
     end
     segIdx = randperm(numSeg,tarObjColocWithRef);
     
     for i = 1:tarObjColocWithRef
         
         segCoordEnd = sum(coordsPerSeg(1:segIdx(i)));
         if segIdx(i) == 1
             segCoordStart = 1;
         else
            segCoordStart = sum(coordsPerSeg(1:segIdx(i)-1)) + 1;
         end
         
        tarColocWithRef(i,:) = datasample(segBoundaries(segCoordStart:segCoordEnd,:),1,'Replace',false);
    
     end
    tarColocWithRef = tarColocWithRef + (colocDistRange(2)-colocDistRange(1))*rand(tarObjColocWithRef,2) + colocDistRange(1);
          
    condColocwithRefImgIdx = sub2ind(size(cellMask),round(tarColocWithRef(:,1)),round(tarColocWithRef(:,2)));
    
     tarWrefImg(condColocwithRefImgIdx) = 1;
 
 end
 if tarBlobOrPt
     tarNotwRefImg = imdilate(tarNotwRefImg,strel('line',lineLength,randi([0 180],1)));
     tarNotwRefImg = logical(imdilate(tarNotwRefImg,strel('disk',tarDiskSize)));
     
     
     tarWrefImg = imdilate(tarWrefImg,strel('line',lineLength,randi([0 180],1)));
     tarWrefImg = logical(imdilate(tarWrefImg,strel('disk',tarDiskSize)));

 end
 
 tarImg = logical(tarNotwRefImg + tarWrefImg);
%  figure, imshowpair(tarImg,segmentation);
%  figure, imshowpair(tarNotwRefImg,segmentation)
%  figure, imshowpair(tarWrefImg,segmentation)
  
%%Run colocalization
[pTR, nullTR] = colocalMeasureBlob2Blob(segmentation,tarImg,cellMask,distThresh,numRandomizations);

end