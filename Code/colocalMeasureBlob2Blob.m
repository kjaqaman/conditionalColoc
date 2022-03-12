function [cT, cNull] = colocalMeasureBlob2Blob(refSeg,tarSeg,cellMask,threshold,numRandomizations)
%COLOCALMEASUREBLOB2BLOB measures colocalization between two sets of segmented objects 
%
%SYNOPSIS: [cT, cNull] = colocalMeasureBlob2Blob(refSeg,tarSeg,cellMask,threshold,numRandomizations)
%Function uses nearest neighbors distances and takes the average minimum 
%distance from target to the reference to measure colocalization. In the
%case target objects are punctate, function will measure colocalization as
%described in function colocalMeasurePt2Blob.m.Punctate objects can be
%input by creating a binary images where ones represent the position of
%each detected object.
%
%INPUT
%   refSeg: binary image of segmented object to be used as the reference
%
%   tarSeg: binary image of segmented object to be used as target
%
%   cellMask: binary mask of ROI (must be all ones if no real mask is input)
%             Needed to determine analyis area
%
%OPTIONAL
%   threshold: distance threshold for target-reference colocalization 
%             Default: 3
%
%   numRandomizations: number of randomizations done to assess the
%                      significance of target-reference colocalization
%                      Default: 1
%
%OUTPUT
%   cT: probability of target colocalizing with reference
%
%   cNull: probability of coincidental target-reference colocalization
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
%% Parse input

%Check reference image is input properly
if isempty(refSeg) || (~islogical(refSeg) && ~isa(refSeg,'double'))
    error('A binary image containing a segmentation must be input for the reference channel')
end

%Check targte image is input properly
if isempty(tarSeg) || (~islogical(tarSeg) && ~isa(tarSeg,'double'))
    error('A binary image containing a segmentation must be input for the target channel')
end

%make sure a mask is input and in the correct format
if isempty(cellMask) || (~islogical(cellMask) && ~isa(cellMask,'double'))
    error('A mask defining the are for analysis must be input')
end

%check threshold
if nargin < 3 || isempty(threshold)
    threshold = 3;
end

%check number of randomizations
if nargin < 4 || isempty(numRandomizations) 
    numRandomizations = 1;
end

%% Analysis

%get object's coordinates
tarSegInfo = regionprops(logical(refSeg),'PixelList','Centroid');

refSegInfo = struct2cell(tarSegInfo)';

tarSegInfo = regionprops(logical(tarSeg),'PixelList','Centroid');
nTarObjects = length(tarSegInfo);

tarSegInfo = struct2cell(tarSegInfo)';

% cT calculation
nnDist = NaN(nTarObjects,1);
refCoords = vertcat(refSegInfo{:,2});
parfor tarBlob = 1:nTarObjects
    
   [~,d] = knnsearch(refCoords,tarSegInfo{tarBlob,2});
   nnDist(tarBlob) = mean(d);    
end

test = find(nnDist <= threshold);
cT = length(test)/nTarObjects;

%% Randomization

blobSizes = cellfun(@length,tarSegInfo(:,2));

if any(blobSizes == 1)%do grid when objects are punctate
    [~, RandDist] = knnsearch(refCoords,cellMask,'K',1);
    
    test = find(RandDist<=threshold);
    cNull = length(test)/length(RandDist);

else%do randomization when objects are not punctate

%preaollocate space for random blob coordinates
condRandCoordVec = cell(nTarObjects,1,numRandomizations);
cNullS = NaN(numRandomizations,1);

for repeat = 1:numRandomizations
    
blobRandShift = NaN(nTarObjects,2);
[~, sortIdx] = sort(blobSizes,'descend');

%create the mask to be erode on each iteration(this keeps the input mask untouched)
    erodeMask = padarray(cellMask,[1 1]);
    %creat mask to store each new blob position. This mask will keep trak
    %of all new blob positions
    cellMaskWithblobs = cellMask; 
    
for i = sortIdx'
        %get the number of pixels to erode the mask. Divided by pi to get the radius
        blobsqrtArea = round(sqrt(blobSizes(i)/pi));
        
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
        blobRandShift(i,2) = blobRandCenter(1) - round(tarSegInfo{i,1}(2));
        blobRandShift(i,1) = blobRandCenter(2) - round(tarSegInfo{i,1}(1));
        
        %shift each blob pixel by blobRandShift (image coord)
        condRandCoordVec{i,1,repeat}(:,1) = tarSegInfo{i,2}(:,1) + blobRandShift(i,1);
        condRandCoordVec{i,1,repeat}(:,2) = tarSegInfo{i,2}(:,2) + blobRandShift(i,2);
        
        %make sure all pixels are with the mask and image limits
        condRandCoordVec{i,1,repeat} = max(condRandCoordVec{i,1,repeat},ones(size(condRandCoordVec{i,1,repeat}))); 
        
        condRandCoordVec{i,1,repeat}(:,1) = min(condRandCoordVec{i,1,repeat}(:,1),size(cellMask,2).*ones(size(condRandCoordVec{i,1,repeat}(:,1))));
        condRandCoordVec{i,1,repeat}(:,2) = min(condRandCoordVec{i,1,repeat}(:,2),size(cellMask,1).*ones(size(condRandCoordVec{i,1,repeat}(:,2))));
        
        %get index of shifted blob coordinates
        shiftBlobCoord = sub2ind(size(cellMask),condRandCoordVec{i,1,repeat}(:,2),condRandCoordVec{i,1,repeat}(:,1));
        %make zero the are where new blob is
        cellMaskWithblobs(shiftBlobCoord) = 0;
        erodeMask = cellMaskWithblobs;
end

%cNull calculation
nnRandDist = NaN(nTarObjects,1);
parfor tarBlob = 1:nTarObjects
    
   [~,d] = knnsearch(refCoords,condRandCoordVec{tarBlob,1,repeat});
   nnRandDist(tarBlob) = mean(d);    
end

randTest = find(nnRandDist <= threshold);
cNullS(repeat) = length(randTest)/nTarObjects;

end
%take mean over all repeats 
cNull = mean(cNullS,'omitnan');
end

end