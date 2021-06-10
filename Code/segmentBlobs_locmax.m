
function [maskBlobs,labels,maskBlobsVis] = segmentBlobs_locmax(image,...
    thresholdMethod,methodValue,filterNoise,filterBackground,minSize,locMax,plotRes,mask) 
%segmentBlobs_locmax segments blobs in 2Detection_IdxD images via various 
%thresholding methods and use local maxima to catch the dim objects
%
%SYNOPSIS [maskBlobs, labels, maskBlobVis] = blobSegmentThreshold(image,...
%  threshold,methodValue,filterNoise,filterBackground,minSize,plotRes,mask)
%
%INPUT  image     : 2D image to be segmented.
%
%       thresholdMethod: 
%                   'otsu'           for Otsu.
%                   'rosin'          for Rosin.
%                   'minmax'         for first minimum after first maximum.
%                   'prct'           to use a certain percentile.
%                   'weightedThresh' weighted combination of Otsu and Rosin
%                   'user'           to use a threshold input by user.
%                   Optional. Default: 'otsu'.
%
%       methodValue: Needed only if thresholdMethod = 'prct' or 'user'.
%                    If 'prct', then this is the percentile to use.
%                    Optional. Default: 90.
%                    If 'user', then this is the threshold value to use.
%                    Optional. Default: 0.9.
%
%       filterNoise: Either 0 to not filter noise or filter sigma > 0 to
%                    filter noise.
%                    Optional. Default: 1.
%
%       filterBackground: Either 0 to not filter background or filter sigma
%                         > 0 to filter background.
%                         Optional. Default: 10.
%
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%
%       locMax    : 1 find local maxima after segmentation. 0 otherwise
%                   Optional. Default 0 
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%
%       mask      : Binary mask. Optional. If not provided, the whole image
%                   domain is segmented.
%
%OUTPUT maskBlobs : mask including segmentations and local maxima.  
%     
%
%      Labels  : Labels of objects in maskBlobs.
%   
%Khuloud Jaqaman October 2012
%Jesus Vega-Lugo July 2017
%Updated by JVL, August 2019. Added maximum filter and weightedThresh
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

%check number of input arguments
if nargin < 1
    disp('Please enter at least image to be segmented');
    return
end

%thresholding method
if nargin < 2 || isempty(thresholdMethod)
    thresholdMethod = 'otsu';
end

%method value
if nargin < 3 || isempty(methodValue)
    switch thresholdMethod
        case 'prct'
            methodValue = 90;
        case 'user'
            methodValue = 0.9;
        otherwise
            methodValue = [];
    end
end

%noise filtering
if nargin < 4 || isempty(filterNoise)
    filterNoise = 1;
end

%background filtering
if nargin < 5 || isempty(filterBackground)
    filterBackground = 10;
end

%minimum blob size
if nargin < 6 || isempty(minSize)
    minSize = 20;
end

%do local maxima
if nargin < 7 || isempty(locMax)
    locMax = 0;
end

%plot results
if nargin < 8 || isempty(plotRes)
    plotRes = 0;
end

%mask
if nargin < 9 || isempty(mask)
    mask = ones(size(image));
end

if ~logical(mask)
    error('Mask must be a logical image.');
end
    
%% Segmentation

%make sure that image is in double format
image = double(image);
mask = logical(mask);
%remove noise by filtering image with a narrow Gaussian
if filterNoise > 0
    imageFiltered = filterGauss2D(image,filterNoise);
else
    imageFiltered = image.*mask;
end 

%estimate background by filtering image with a wide Gaussian
if filterBackground > 0
    imageBackground = filterGauss2D(image,filterBackground);
else
    imageBackground = zeros(size(image));
end

%calculate noise-filtered and background-subtracted image
imageFilteredMinusBackground = imageFiltered - imageBackground;

imageFilteredMinusBackground = imageFilteredMinusBackground .* mask;
%enhance features by performing a maximum filter
[sizeX,sizeY] = size(imageFilteredMinusBackground);
imageDilated = imageFilteredMinusBackground;
imageTmp(:,:,1) = imageDilated;
imageTmp(:,:,2) = [zeros(1,sizeY); imageDilated(2:end,:)];
imageTmp(:,:,3) = [imageDilated(1:end-1,:); zeros(1,sizeY)];
imageTmp(:,:,4) = [zeros(sizeX,1) imageDilated(:,2:end)];
imageTmp(:,:,5) = [imageDilated(:,1:end-1) zeros(sizeX,1)];
imageTmp(:,:,6) = [zeros(1,sizeY); [zeros(sizeX-1,1) imageDilated(2:end,2:end)]];
imageTmp(:,:,7) = [zeros(1,sizeY); [imageDilated(2:end,1:end-1) zeros(sizeX-1,1)]];
imageTmp(:,:,8) = [[zeros(sizeX-1,1) imageDilated(1:end-1,2:end)]; zeros(1,sizeY)];
imageTmp(:,:,9) = [[imageDilated(1:end-1,1:end-1) zeros(sizeX-1,1)]; zeros(1,sizeY)];
imageDilated = max(imageTmp,[],3);

%find nonzero values (due to masking)
nzInd = find(imageDilated);

%get minumum and maximum pixel values in image
minSignal = min(imageDilated(nzInd));
maxSignal = max(imageDilated(nzInd));

%normalize nonzero value between 0 and 1
imageDilatedNorm = zeros(size(imageDilated));
imageDilatedNorm(nzInd) = (imageDilated(nzInd) - minSignal) / (maxSignal - minSignal);

 
%estimate the intensity level to use for thresholding the image
switch thresholdMethod
    case 'otsu'
        try
            level = graythresh(imageDilatedNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (otsu, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'rosin'
        try
            [dummy,level] = cutFirstHistMode(imageDilatedNorm(nzInd),0); %#ok<ASGLU>
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (rosin, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'minmax'
        try
            level = thresholdFluorescenceImage(imageDilatedNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (minmax, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'prct'
        try
            level = prctile(imageDilatedNorm(nzInd),methodValue);
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (prct, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'user'
        try
            level = methodValue;
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (user, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'weightedThresh'
        try 
            level1 = graythresh(imageDilatedNorm(nzInd)); %Otsu
            [~, level2] = cutFirstHistMode(imageDilatedNorm(nzInd),0); %Rosin
            level = 0.33333*level1 + 0.66667*level2;
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (weightedThresh,'...
                'noise ' num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
end

imageThresholded = imbinarize(imageDilatedNorm,level); 

%fill holes in thresholded image to make continuous blobs
imageThresholdedFilled = imfill(imageThresholded,'holes');

% go over blobs and remove those with a size smaller that minSize
labels = bwlabel(imageThresholdedFilled);
labels = labels.*mask;
%Take the objects with and area bigger than minSize
%JVL March 2019:removed logical from following line
stats = regionprops(labels, 'Area'); 
idx = find([stats.Area] > minSize);

%mask with blobs
maskBlobs = ismember(labels, idx);

SE = strel('disk',1);
maskBlobsVis = imdilate(maskBlobs,SE); 

%% Local maxima
%get local maxima

if locMax
locmaxvalues = locmax2d(imageDilatedNorm,3,[]);
%locmaxvalues = locmax2d((image-min(image(:)))/(max(image(:))-min(image(:))),3,[]);

% get position of local maxima
locmaxidx = find(locmaxvalues);

%get local maxima inside blobs
locmaxInBlobs = ismember(locmaxidx,find(maskBlobs));

%remove local maxima that are inside a blob
% this eliminates bright spot that are already included in maskBlobs
locmaxidx(locmaxInBlobs) = [];

%get values inside mask
locmaxInMask = ismember(locmaxidx,find(mask));

%get actual value of the local maxima
locmaxValuesInMask = locmaxvalues(locmaxidx(locmaxInMask));


%take the mean of the local maxima excluding outliers
[locmaxMean, locmaxStd] = robustMean(locmaxValuesInMask,[],3);

% %take values above the third std 
threshValue = locmaxMean+3*locmaxStd;
 
% %threshold local maxima
 thresh = imbinarize(locmaxvalues,threshValue);
 
%include blobs and local maxima in the same mask
%threshidx= find(thresh);%find(locmaxvalues>=thresh);
maskBlobs(thresh>=thresh) = 1;
maskBlobs = maskBlobs.*mask; 
 SE = strel('disk',1);
 maskBlobsVis = imdilate(maskBlobs,SE); 
 %get labels of mask with segementation and local maxima
 labels = bwlabel(maskBlobs);
end

%% Plotting
 if plotRes
    
    imageScaled = (image - prctile(image(:),0.01)) / (prctile(image(:),99.99) - prctile(image(:),0.01));
    imageScaled(imageScaled<0) = 0;
    imageScaled(imageScaled>1) = 1;
    
    %get the blob edges from the final blob mask
    SE = strel('square',3);
    edgesBlobs = imdilate(maskBlobs,SE) - maskBlobs;
    
    %give the edge pixels a value of zero in the original image
    imageScaled(edgesBlobs==1) = 0;
    
    
    %construct a 3-layered image to show blob edges on top of
    %original image
    image3Color = repmat(imageScaled,[1 1 3]);
    image3Color(:,:,1) = image3Color(:,:,1) + edgesBlobs;
    
         
    figure('Name',['segmentation_' thresholdMethod '_noise' num2str(filterNoise)...
       '_background' num2str(filterBackground) '_min=' num2str(minSize)]);
        imshow(image3Color,[]);
       
   
 

end
