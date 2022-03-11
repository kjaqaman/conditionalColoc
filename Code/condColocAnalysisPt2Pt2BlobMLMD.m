function condColocAnalysisPt2Pt2BlobMLMD(MLMD,detectedChannels, blobChannel,varargin)
%CONDCOLOCANALYSISPT2PT2BlobMLMD run conditional colocalization analysis for 2 punctate objects and 1 non-punctate object. 
%
%SYNOPSIS condColocAnalysisPt2Pt2BlobMLMD(MLMD,detectedChannels,blobChannel,varargin)
%
% Function will run conditional colocalization analysis for two punctate
% objects and one non-punctate object. Punctate objects will permutated to
% be target, reference, and condition; non-punctate object will be
% permutated to be reference and condition.
%
%INPUT 
%   Required
%       MLMD:               MovieList or MovieData object for movie(s) to be analyzed
%
%       detectedChannels:   Vector containing the indices of the two
%                           punctate channels.  e.g. [1 2]
%                              
%
%       blobChannel:        non-punctate object channel index
%                           
%   Optional input arguments as name-value pairs
%       detectionProcessID: Cell array containing desired detection process 
%                           to be used for each channel with punctate objects.
%                           E.g.{DetectProcessForChA DetectProcessForChB}
%                           Options: 'SubRes'and 'PointSource'
%                           NOTE: SubRes refers to "Gaussian mixture-model
%                           fitting" in u-Track.
%                           Default: 'SubRes'(for both channels)
%                           i.e. {'SubRes' 'SubRes'}
%
%       colocDistThresh:    Vector containing colocalization distance threshold 
%                           for pairs of objects. 
%                           colocDistThresh(1): threshold between 
%                                               punctate objects
%                           colocDistThresh(2): threshold between
%                                               the first punctate channel
%                                               (detectedChannels(1)) and
%                                               blobChannel.
%                           colocDistThresh(3): threshold between
%                                               the second punctate channel
%                                               (detectedChannels(2)) and
%                                               blobChannel.
%                           E.g. Let detectedChannels = [2 3] and
%                           blobChannel = 1. If colocDistThresh = [5 1 2]
%                           then the distance threshold between 2 and 3 =
%                           5, threshold between 2 and 1 = 1, and threshold
%                           between 3 and 1 = 2.
%                           Default: [3 3 3] (pixels)
%
%   numRandomizations: number of randomizations to be used for calculating
%                      randC value. (see vega-Lugo et al. for randC defintion)
%                      Default: 100
%
%       alphaValue:         Alpha value for comparing real data to random
%                           (see colocalMeasurePt2Pt)
%                           Default: 0.05
%                           NOTE: for conditional colocalization analysis
%                           as described in Vega-Lugo et al. 2022 this
%                           parameter is not relevant.
%
%OUPUT Output is saved in directory 
%      For detailed list of output parameters see condColocPt2Pt2Blob.m documentation
%
%Jesus Vega-Lugo June, 2019
%
%Jesus Vega-Lugo Dec., 2020 modified to use different distance thresholds
%for different pair of molecules. Before same threshold was used for all
%pairs of molecules.
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


ip = inputParser;

ip.CaseSensitive = false;
ip.FunctionName  = 'condColocAnalysisPt2Pt2BlobMLMD';

%Add required parameter
addRequired(ip, 'MLMD', @(x) isa(x,'MovieData') || isa(x,'MovieList'))

addRequired(ip,'detectedChannels',@isvector)

addRequired(ip,'blobChannel', @isscalar)

%Add optional parameters
addParameter(ip,'detectionProcessID',{'SubRes','SubRes'}, @iscell)
addParameter(ip,'colocDistThresh', [3 3 3], @isvector)
addParameter(ip,'numRandomizations',100,@isscalar)
addParameter(ip,'alphaValue', 0.05, @isscalar)

parse(ip,MLMD,detectedChannels,blobChannel,varargin{:})


%Check there is no conflict between detected and blob channel input
if any(detectedChannels == blobChannel)
    error('blobChannel must be different from both indices in detectedChannels')
end

%% Analysis

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

%go over all movies and run Pt2Pt colocalization
for i = 1 : numMovies
    
    %get movieData of current movie
    if listFlag == 1
        MD = MovieData.load(ML.movieDataFile_{i});
    end
    
    iProcCondColoc = MD.getProcessIndex('CondColocPt2Pt2BlobProcess',Inf,0);
    
    if isempty(iProcCondColoc)
       iProcCondColoc = numel(MD.processes_) + 1;
       MD.addProcess(CondColocPt2Pt2BlobProcess(MD));
    end
        
    
    %define colocalization parameters
    
    p = MD.getProcess(iProcCondColoc).getParameters();
   
    p.DetectedChannels = detectedChannels;
    p.ChannelBlob = blobChannel;
    p.DetectionProcessID = ip.Results.detectionProcessID;
        
    p.ColocDistThresh = ip.Results.colocDistThresh;
    p.NumRandomizations = ip.Results.numRandomizations;
    p.AlphaValue = ip.Results.alphaValue;
    
    MD.getProcess(iProcCondColoc).setParameters(p);
    MD.getProcess(iProcCondColoc).run
    
   
end

end
