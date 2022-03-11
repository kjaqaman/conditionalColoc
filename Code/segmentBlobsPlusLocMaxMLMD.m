function segmentBlobsPlusLocMaxMLMD(MLMD,blobChannel,varargin)
%SEGMENTBLOBSPLUSLOCMAXMLMD runs segmentation process SegmentBlobsPlusLocMaxSimpleProcess
%
%Synopsis segmentBlobsPlusLocMaxMLMD(MLMD,varargin)
%
%INPUT
%Mandatory:
%   MLMD: MovieData or MovieList to be segmented
%       
%       blobChannel: scalar or vector containing channel index(ces) for blobs
%                   
%Optional:
%       thresholdMethod: 
%                   'otsu'           for Otsu.
%                   'rosin'          for Rosin.
%                   'minmax'         for first minimum after first maximum.
%                   'prct'           to use a certain percentile.
%                   'weightedThresh' weighted combination of Otsu and Rosin
%                   'user'           to use a threshold input by user.
%                   Optional. Default: 'prct'.
%
%       methodValue: Needed only if thresholdMethod = 'prct' or 'user'.
%                    If 'prct', then this is the percentile to use.
%                    If 'user', then this is the threshold value to use.
%                    Optional. Default: 90
%
%       filterNoise: Either 0 to not filter noise or filter sigma > 0 to
%                    filter noise.
%                    Optional. Default: 1.
%
%       filterBackground: Either 0 to not filter background or filter sigma
%                         > 0 to filter background.
%                         Optional. Default: 50.
%
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%
%       locMax    : 1 find local maxima after segmentation. 0 otherwise
%                   Default: 0.
%
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Default: 0.
%
%OUTPUT maskBlobs : mask including segmentations and local maxima.
%   Results store in directory
%
%Jesus Vega-Lugo November 2019
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
%% Inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.FunctionName = 'segmentBlobsPlusLocMaxMLMD';

%require moviedata or movielist
addRequired(ip,'MLMD',@(x) isa(x,'MovieData') || isa(x,'MovieList'))

addRequired(ip,'blobChannel',@(x) isscalar(x) || isvector(x))

%optional parameters
addParameter(ip,'thresholdMethod','prct',@(x) ischar(x) || iscell(x))
addParameter(ip,'methodValue',90,@(x) isscalar(x) || isvector(x))
addParameter(ip,'filterNoise',1,@(x) isscalar(x) || isvector(x))
addParameter(ip,'filterBackground',50,@(x) isscalar(x) || isvector(x))
addParameter(ip,'minSize',20,@(x) isscalar(x) || isvector(x))
addParameter(ip,'locMax',0,@(x) isscalar(x) || isvector(x))
addParameter(ip,'plotRes',0,@(x) isscalar(x) || isvector(x))

parse(ip,MLMD,blobChannel,varargin{:})

blobChannel = ip.Results.blobChannel;
thresholdMethod = ip.Results.thresholdMethod;
methodValue = ip.Results.methodValue;
filterNoise = ip.Results.filterNoise;
filterBackground = ip.Results.filterBackground;
minSize = ip.Results.minSize;
locMax = ip.Results.locMax;
plotRes = ip.Results.plotRes;

if isvector(blobChannel)
    nChan = length(blobChannel);
    %when multiple channels are to be segmented check that each paramater
    %is input as a vector 
    checkParams = [~iscell(thresholdMethod), isscalar(methodValue), isscalar(filterNoise),...
        isscalar(filterBackground), isscalar(minSize),isscalar(locMax), isscalar(plotRes)];
    
    %if specified as scalar, they will be replicated
    for i = 1:7
        if checkParams(i)
           switch i
               case 1
                   t = thresholdMethod;
                   thresholdMethod  = cell(1,nChan);
                   thresholdMethod(cellfun(@isempty,thresholdMethod)) = {t};
               
               case 2
                   t = methodValue;
                   methodValue = NaN(1,nChan);
                   for ii = 1:nChan
                       switch thresholdMethod{ii}
                           case {'prct' 'user'}
                               methodValue(1,ii) = t;
                       end
                   end
                   
               case 3
                   filterNoise(1,1:nChan) = filterNoise;
                   
               case 4
                   filterBackground(1,1:nChan) = filterBackground;
                   
               case 5
                   minSize(1,1:nChan) = minSize;
                   
               case 6
                   locMax(1,1:nChan) = locMax;
                   
               case 7 
                   plotRes(1,1:nChan) = plotRes;
           end
        end
    end%for 
end%if chIndex is vector

%% Run movies

%determine if input is a MovieList or MovieData object
if isa(MLMD,'MovieList') %if it is a movieList
    
    listFlag = 1;
    
    %rename to ML
    ML = MLMD;
    clear MLMD
    
    %get number of movies
    numMovies = length(ML.movieDataFile_);
    
else %if it is a movieData
    
    listFlag = 0;
    
    %rename to MD
    MD = MLMD;
    clear MLMD
    
    %assign number of movies to 1
    numMovies = 1;
    
end

%go over all movies
for iM = 1 : numMovies
    %load movie from movielist
    if listFlag == 1
            MD = MovieData.load(ML.movieDataFile_{iM});
    end
    
%% Run segmentation

    %look for existibg process
    iSeg = MD.getProcessIndex('SegmentBlobsPlusLocMaxSimpleProcess',1,0);
    
    %if it does not exis, creste it
    if isempty(iSeg)
        iSeg = length(MD.processes_) + 1;
        MD.addProcess(SegmentBlobsPlusLocMaxSimpleProcess(MD))
    end
    
    %setparameters
    pSeg = MD.getProcess(iSeg).getParameters();
    
    if ~iscell(pSeg.detectionParam.thresholdMethod)
        pSeg.detectionParam.thresholdMethod = cell(1,nChan);
    end
    
    pSeg.ChannelIndex = blobChannel;
    pSeg.detectionParam.thresholdMethod(blobChannel) = thresholdMethod;
    pSeg.detectionParam.methodValue(blobChannel) = methodValue;
    pSeg.detectionParam.filterNoise(blobChannel) = filterNoise;
    pSeg.detectionParam.filterBackground(blobChannel) = filterBackground;
    pSeg.detectionParam.minSize(blobChannel) = minSize;
    pSeg.detectionParam.locMax(blobChannel) = locMax;
    pSeg.plotRes(blobChannel) = plotRes;
    
    %run process
    MD.getProcess(iSeg).setParameters(pSeg)
    MD.getProcess(iSeg).run
end

end