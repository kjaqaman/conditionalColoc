function colocMeasure = simulateCondColocGeneral(numRefObj,numTarObj,segmentation,...
                cellMask,refColocWithCondFrac,tarColocWithCondFrac,...
                tarColocWithRefAtCondFrac,tarColocWithRefNotAtCondFrac,varargin)
%SIMULATEPT2PT2BLOBCONDCOLOC simulate object positions to be used for general conditional coloalization
%
% SYNOPSIS colocMeasure, numOfObjects]  = simulatePt2Pt2BlobCondColoc(numRefDetect,...
%                numTarDetect,sementation,cellMask,refColocWithCondFrac,tarColocWithCondFrac,...
%                tarColocWithRefAtCondFrac,tarColocWithRefNotAtCondFrac,varargin)
%
% Function sets up simulated positions and runs colocalization general 
% conditional colocalization analysis. When creating non-punctate target or 
% reference objects, function will first place the center of the the object
% then it will dilate the center with a line of specified length (see
% lineLength) and further dilated with a disk of specified radius (see 
% tarDiskSize and refDiskSize)
%
%INPUT
% Required
%   numRefObj: number of reference objects
%
%   numTarObj: number of traget objects
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
%
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
%
%   numTarRefRand: number of target randomizations to assess conditional
%                  colocalization between target and reference
%                  Default 1:
%
%   numCondRand: number of condition randomizations to calculate randC
%                (see vega-Lugo et al. 2022 for randC defintion)
%                Default: 1 
%
%   tarBlobOrPt: 1 to create non-punctate traget object, 0 for creating
%                punctate target object. Default: 1
%
%
%   refBlobOrPt: 1 to create non-punctate reference object, 0 for creating
%                punctate reference object. Default: 1
%
% Below parameters are only needed if tarBlobOrPt and/or refBlobOrPt = 1
%
%   lineLength: length of line to be used to dilate the center of
%               non-punctate target and reference objects. Default: 7
%
%   tarDiskSize: radius for disk to dilate initial line dilation of target
%                objects. Default: 5
%
%   refDiskSize: radius for disk to dilate initial line dilation of
%                reference objects. Default: 3
%
%OUTPUT
%   colocalMeasure: see output documentation on condColcoPt2Pt2Blob.m
%
%Jesus Vega-Lugo Dec. 2021
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
addParameter(ip,'numTarRefRand',1,@isscalar)
addParameter(ip,'numCondRand',1,@isscalar)
addParameter(ip,'tarBlobOrPt',1,@(x) x == 0 || x == 1)
addParameter(ip,'refBlobOrPt',1,@(x) x == 0 || x == 1)
addParameter(ip,'refDiskSize',5,@isscalar)
addParameter(ip,'tarDiskSize',3,@isscalar)
addParameter(ip,'lineLength',7,@isscalar)

parse(ip,numRefObj,numTarObj,segmentation,cellMask,refColocWithCondFrac,...
        tarColocWithCondFrac,tarColocWithRefAtCondFrac,...
        tarColocWithRefNotAtCondFrac,varargin{:})


numRefObj = ip.Results.numRefDetect;
numTarObj = ip.Results.numTarDetect;
segmentation = ip.Results.segmentation;
cellMask = ip.Results.cellMask;
refColocWithCondFrac = ip.Results.refColocWithCondFrac;
tarColocWithCondFrac = ip.Results.tarColocWithCondFrac;
tarColocWithRefAtCondFrac = ip.Results.tarColocWithRefAtCondFrac;
tarColocWithRefNotAtCondFrac = ip.Results.tarColocWithRefNotAtCondFrac;
colocDistRange = ip.Results.colocDistRange;
colocDistThresh = ip.Results.colocDistThresh;
numTarRefRand = ip.Results.numTarRefRand;
numCondRand = ip.Results.numCondRand;
tarBlobOrPt = ip.Results.tarBlobOrPt;
refBlobOrPt = ip.Results.refBlobOrPt;
refDiskSize = ip.Results.refDiskSize;
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
                        
%%    %%%%%%%%%%%%%%%%%%%%%%Reference%%%%%%%%%%%%%%%%%%%%%%%
imgSize = size(cellMask);

refWcondImg = zeros(imgSize);
refNotwCondImg = zeros(imgSize);

%number of colocalized ref molecules
 refDetectColocWithCond = round(numRefObj*refColocWithCondFrac);

 numRefNotCond = numRefObj-refDetectColocWithCond;
 
 refNotColocWithCond = NaN(numRefNotCond,2);
 
 refPtMask = erodedCellMask;
 
 %place ref not colocalizing with cond
 for i = 1:numRefNotCond
    refNotColocWithCond(i,:) = datasample(maskWithoutBlobsCoords,1,'Replace',false);
    maskWithoutBlobsCoords = [];
    refPtMask(refNotColocWithCond(i,1),refNotColocWithCond(i,2)) = 0;
    refMask = imerode(refPtMask,strel('disk',refDiskSize+5));
    refMaskWithoutBlobs = refMask.*erodedMask;
    [maskWithoutBlobsCoords(:,1), maskWithoutBlobsCoords(:,2)] = find(refMaskWithoutBlobs);
 end   

 refNotColocWithCondIdx = sub2ind(size(cellMask),refNotColocWithCond(:,1),refNotColocWithCond(:,2));
 
 refNotwCondImg(refNotColocWithCondIdx) = 1;

%place ref colocalizing with condition
 if refDetectColocWithCond
     refColocWithCond = NaN(refDetectColocWithCond,2);
     numSeg = length(coordsPerSeg);    
     segIdx = randperm(numSeg,refDetectColocWithCond);
     for i = 1:refDetectColocWithCond
         
         segCoordEnd = sum(coordsPerSeg(1:segIdx(i)));
         if segIdx(i) == 1
             segCoordStart = 1;
         else
            segCoordStart = sum(coordsPerSeg(1:segIdx(i)-1)) + 1;
         end
         
        refColocWithCond(i,:) = datasample(segBoundaries(segCoordStart:segCoordEnd,:),1,'Replace',false);
    
     end
    refColocWithCond = refColocWithCond + (0.2*rand(refDetectColocWithCond,2)+0.2);
    condColocwithRefImgIdx = sub2ind(size(cellMask),round(refColocWithCond(:,1)),round(refColocWithCond(:,2)));
    
     refWcondImg(condColocwithRefImgIdx) = 1;
     
 else
    refColocWithCondPositions = [];
 end

 %dilate ref centers when specified as non-punctate
 if refBlobOrPt
     refNotwCondImg = imdilate(refNotwCondImg,strel('line',lineLength,randi([0 180],1)));
     refNotwCondImg = logical(imdilate(refNotwCondImg,strel('disk',refDiskSize)));
     
     refNotColocWithCondPositions = regionprops(logical(refNotwCondImg),'PixelList');
     refNotColocWithCondPositions = cell2mat(struct2cell(refNotColocWithCondPositions)');
     
     
     refWcondImg = imdilate(refWcondImg,strel('line',lineLength,randi([0 180],1)));
     refWcondImg = logical(imdilate(refWcondImg,strel('disk',refDiskSize)));
     
     refColocWithCondPositions = regionprops(logical(refWcondImg),'PixelList');
     refColocWithCondPositions = cell2mat(struct2cell(refColocWithCondPositions)');
 end
 
 refImg = logical(refNotwCondImg + refWcondImg);
%   figure, imshowpair(refImg,segmentation);
%   figure, imshowpair(refNotwCondImg,segmentation)
%   figure, imshowpair(refWcondImg,segmentation)
  

refColocWithCondPositions = refColocWithCondPositions(:,[2 1]);
refNotColocWithCondPositions = refNotColocWithCondPositions(:,[2 1]);
  
allRefPositions = vertcat(refColocWithCondPositions,refNotColocWithCondPositions);

%%    %%%%%%%%%%%%%%%%%%%%%%%Target%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number tar molecules colocalizing with condition +- ref
tarDetectColocWithCond = round(numTarObj*tarColocWithCondFrac);
tarDetectColocWithRefAtCond = round(tarDetectColocWithCond*tarColocWithRefAtCondFrac);

%%%%%%%%%%%%%%%%%%%%%% tar + cond + ref%%%%%%%%%%%%%%%%%%%%
if tarDetectColocWithRefAtCond %tar colocalizing with cond and ref
    
    numRefColocWithCond = length(refColocWithCond);
    if numRefColocWithCond >= tarDetectColocWithRefAtCond
        
        
        %using ref centroid to place tar
        
        tarDetectColocWithRefAtCondPositions = refColocWithCond(randperm(numRefColocWithCond,tarDetectColocWithRefAtCond),:)...
            + ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefAtCond,2) + colocDistRange(1));
        
    else
        
        %using ref centroid to palce tar
        repRefColocWithCondPositions = refColocWithCond;
        while  length(repRefColocWithCondPositions) < tarDetectColocWithRefAtCond
            repRefColocWithCondPositions = repmat(repRefColocWithCondPositions,2,1);
        end
        numrep = length(repRefColocWithCondPositions);
        tarDetectColocWithRefAtCondPositions = refColocWithCond(randperm(numrep,tarDetectColocWithRefAtCond),:)...
            + ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefAtCond,2) + colocDistRange(1));
    end
else
    tarDetectColocWithRefAtCondPositions = [];
end


tarDetectColocWithRefAtCondIdx = sub2ind(size(cellMask),...
    round(tarDetectColocWithRefAtCondPositions(:,1)),round(tarDetectColocWithRefAtCondPositions(:,2)));

 tarDetectColocWithRefAtCondImg = zeros(imgSize);
 tarDetectColocWithRefAtCondImg(tarDetectColocWithRefAtCondIdx) = 1;
 
 
 if tarBlobOrPt
    tarDetectColocWithRefAtCondImg = imdilate(tarDetectColocWithRefAtCondImg,strel('line',lineLength,randi([0 180],1)));
    tarDetectColocWithRefAtCondImg = logical(imdilate(tarDetectColocWithRefAtCondImg,strel('disk',tarDiskSize)));
  
 end
 
%   tarImg = logical(tarDetectColocWithRefAtCondImg);
%   merge3 = cat(3,double(tarImg),double(refWcondImg),double(segmentation));
%   figure, imshow(merge3);
    

%%%%%%%%%%%%%%%%%%%%%% tar + cond - ref%%%%%%%%%%%%%%%%%%%%

tarDetectColocWithCondNotWithRef = tarDetectColocWithCond - tarDetectColocWithRefAtCond;

segBoundariesNoRefColoc = segBoundaries;

%make sure there are enough conditions not colocalizing with reference
if length(segBoundariesNoRefColoc) >= tarDetectColocWithCondNotWithRef && ~isempty(segBoundariesNoRefColoc)
    if numSeg > tarDetectColocWithCondNotWithRef
       segIdx = randperm(numSeg,tarDetectColocWithCondNotWithRef);
    else
        moreNeeded = tarDetectColocWithCondNotWithRef - numSeg;
        segIdx = horcat(1:numSeg,randperm(numSeg,moreNeeded));
    end
    
    condNotColocWithRef = NaN(tarDetectColocWithCondNotWithRef,2);
    
    for i = 1:tarDetectColocWithCondNotWithRef
        segCoordEnd = sum(coordsPerSeg(1:segIdx(i)));
        if segIdx(i) == 1
            segCoordStart = 1;
        else
            segCoordStart = segCoordEnd - sum(coordsPerSeg(1:segIdx(i)-1)) + 1;
        end
        
        condNotColocWithRef(i,:) = datasample(segBoundariesNoRefColoc(segCoordStart:segCoordEnd,:),1);
    end
    
    condNotColocWithRef = condNotColocWithRef + ((colocDistRange(2)-colocDistRange(1))*...
            rand(tarDetectColocWithCondNotWithRef,2) + colocDistRange(1));
    
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

tarColocWithCondNotWithRefIdx = sub2ind(size(cellMask),...
    round(tarColocWithCondNotWithRefPositions(:,1)),round(tarColocWithCondNotWithRefPositions(:,2)));

 tarColocWithCondNotWithRefImg = zeros(imgSize);
 tarColocWithCondNotWithRefImg(tarColocWithCondNotWithRefIdx) = 1;
 
 if tarBlobOrPt
     tarColocWithCondNotWithRefImg = imdilate(tarColocWithCondNotWithRefImg,strel('line',lineLength,randi([0 180],1)));
    tarColocWithCondNotWithRefImg = logical(imdilate(tarColocWithCondNotWithRefImg,strel('disk',tarDiskSize)));
  
 end
 
%   tarImg = logical(tarColocWithCondNotWithRefImg);
%   merge3 = cat(3,double(tarImg),double(refImg),double(segmentation));
%   figure, imshow(merge3);

tarDetectNotColocWithCond = numTarObj-tarDetectColocWithCond;

tarDetectColocWithRefNotAtCond = round(tarDetectNotColocWithCond*tarColocWithRefNotAtCondFrac);
tarDetectNotColocWithRefOrCond = tarDetectNotColocWithCond - tarDetectColocWithRefNotAtCond;


%%%%%%%%%%%%%%%%%%%%%% tar - cond + ref%%%%%%%%%%%%%%%%%%%%
if tarDetectColocWithRefNotAtCond
    
    numRefNotColocWithCond = length(refNotColocWithCond);
    if numRefNotColocWithCond >= tarDetectColocWithRefNotAtCond
        
        %using ref centroid to place tar
        
        tarColocWithRefNotWithCondPositions = refNotColocWithCond(randperm(numRefNotColocWithCond,tarDetectColocWithRefNotAtCond),:) + ...
            ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefNotAtCond,2) + colocDistRange(1));
    else
        
        repRefNotColocWithCond = refNotColocWithCond;
        while length(repRefNotColocWithCond) < tarDetectColocWithRefNotAtCond
            repRefNotColocWithCond = repmat(repRefNotColocWithCond,2,1);
        end
        numrep = length(repRefNotColocWithCond);
        tarColocWithRefNotWithCondPositions = repRefNotColocWithCond(randperm(numrep,tarDetectColocWithRefNotAtCond),:) + ...
            ((colocDistRange(2)-colocDistRange(1))*rand(tarDetectColocWithRefNotAtCond,2) + colocDistRange(1));
    end
else
    tarColocWithRefNotWithCondPositions = [];
end

tarColocWithRefNotWithCondIdx = sub2ind(size(cellMask),...
    round(tarColocWithRefNotWithCondPositions(:,1)),round(tarColocWithRefNotWithCondPositions(:,2)));

tarColocWithRefNotWithCondImg = zeros(imgSize);
tarColocWithRefNotWithCondImg(tarColocWithRefNotWithCondIdx) = 1;

if tarBlobOrPt
    tarColocWithRefNotWithCondImg = imdilate(tarColocWithRefNotWithCondImg,strel('line',lineLength,randi([0 180],1)));
    tarColocWithRefNotWithCondImg = logical(imdilate(tarColocWithRefNotWithCondImg,strel('disk',tarDiskSize)));   
end
 
%   tarImg = logical(tarColocWithRefNotWithCondImg);
%   merge3 = cat(3,double(tarImg),double(refNotwCondImg),double(segmentation));
%   figure, imshow(merge3);

%%%%%%%%%%%%%%%%tar - cond - ref %%%%%%%%%%%%%%%%%%%

tarPtMask = erodedCellMask;
tarDetectNotColocWithRefOrCondPositions = NaN(tarDetectNotColocWithRefOrCond,2);

for i = 1:tarDetectNotColocWithRefOrCond
   tarDetectNotColocWithRefOrCondPositions(i,:) = datasample(maskWithoutBlobsCoords,1); 
   maskWithoutBlobsCoords = [];
   tarPtMask(tarDetectNotColocWithRefOrCondPositions(i,1),...
       tarDetectNotColocWithRefOrCondPositions(i,2)) = 0;
   tarMask = imerode(tarPtMask,strel('disk',tarDiskSize+1));
   tarMaskWithouBlobs = tarMask.*refMaskWithoutBlobs ;
   [maskWithoutBlobsCoords(:,1), maskWithoutBlobsCoords(:,2)] = find(tarMaskWithouBlobs);
end

tarDetectNotColocWithRefOrCondIdx = sub2ind(size(cellMask),...
    round(tarDetectNotColocWithRefOrCondPositions(:,1)),round(tarDetectNotColocWithRefOrCondPositions(:,2)));

 tarDetectNotColocWithRefOrCondImg = zeros(imgSize);
 tarDetectNotColocWithRefOrCondImg(tarDetectNotColocWithRefOrCondIdx) = 1;
 
 if tarBlobOrPt
     tarDetectNotColocWithRefOrCondImg = imdilate(tarDetectNotColocWithRefOrCondImg,strel('line',lineLength,randi([0 180],1)));
    tarDetectNotColocWithRefOrCondImg = logical(imdilate(tarDetectNotColocWithRefOrCondImg,strel('disk',tarDiskSize)));
   
 end
 
%     tarImg = logical(tarDetectNotColocWithRefOrCondImg);
%     merge3 = cat(3,double(tarImg),double(refImg),double(segmentation));
%     figure, imshow(merge3);
 
   tarImg = logical(tarDetectColocWithRefAtCondImg + tarColocWithCondNotWithRefImg + ...
       tarColocWithRefNotWithCondImg + tarDetectNotColocWithRefOrCondImg);

%% Run conditional colocalization


% merge3 = cat(3,double(tarImg),double(refImg),double(segmentation));
% figure, imshow(merge3)
% bound = bwboundaries(cellMask);
% hold on, plot(bound{1}(:,2),bound{1}(:,1),'r')
% pause(0.5)

colocMeasure = condColocGeneral(refImg,tarImg,segmentation,...
                            cellMask,[colocDistThresh(2) colocDistThresh(3)],colocDistThresh(1),numTarRefRand,numCondRand);

end
