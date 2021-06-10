function segmentBlobsPlusLocMaxMLMD(MLMD,blobChannel,varargin)
%SEGMENTBLOBSPLUSLOCMAXMLMD runs segmentation process SegmentBlobsPlusLocMaxSimpleProcess
%
%Synopsis segmentBlobsPlusLocMaxMLMD(MLMD,varargin)
%
%INPUT
%Mandatory:
%   MLMD: MovieData or MovieList to be segmented
%       
%       blobChannel: channel index for blobs
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

addRequired(ip,'blobChannel',@isscalar)

%optional parameters
addParameter(ip,'thresholdMethod','prct',@ischar)
addParameter(ip,'methodValue',90,@isscalar)
addParameter(ip,'filterNoise',1,@isscalar)
addParameter(ip,'filterBackground',50,@isscalar)
addParameter(ip,'minSize',20,@isscalar)
addParameter(ip,'locMax',0,@isscalar)
addParameter(ip,'plotRes',0,@isscalar)

parse(ip,MLMD,blobChannel,varargin{:})



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

    pSeg.ChannelIndex = ip.Results.blobChannel;
    pSeg.detectionParam.thresholdMethod = ip.Results.thresholdMethod;
    pSeg.detectionParam.methodValue = ip.Results.methodValue;
    pSeg.detectionParam.filterNoise = ip.Results.filterNoise;
    pSeg.detectionParam.filterBackground = ip.Results.filterBackground;
    pSeg.detectionParam.minSize = ip.Results.minSize;
    pSeg.detectionParam.locMax = ip.Results.locMax;
    pSeg.plotRes = ip.Results.plotRes;
    
    %run process
    MD.getProcess(iSeg).setParameters(pSeg)
    MD.getProcess(iSeg).run
end

end