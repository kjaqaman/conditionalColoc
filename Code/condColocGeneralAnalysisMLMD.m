function condColocGeneralAnalysisMLMD(MLMD,ptOrBlob,varargin)
%CONDCOLOCGENERALANALYSIS run colocalization analysis images with three channels with segmented and/or detected objects
%
% Function will run conditional colocalization analysis for three punctate
% or non-punctate objects or a combination. It will make all permutations 
%of target, reference and condition.
% 
%INPUT 
%   Required
%       MLMD:               MovieList or MovieData object for movie(s) to be analyzed
%
%       ptOrBlob:           vector to identify each channel as punctate or
%                           non-punctate. 0 for non-punctate, 1 for
%                           punctate
%
%   Optional input arguments as name-value pairs
%       detectionProcessID: Cell array containing desired detection process 
%                           to be used for each channel with punctate objects.
%                           E.g.{DetectProcessForCh1 DetectProcessForCh2 DetectProcessForCh3}
%                           Options: 'SubRes'and 'PointSource'
%                           NOTE: SubRes refers to "Gaussian mixture-model
%                           fitting" in u-Track.
%                           Default: 'PointSource'(for channel 1)
%                           i.e. {'PointSource', [], []}
%                           NOTE: Non-punctate channel must contain empty
%                           brackets
%
%       colocDistThresh:    Vector containing colocalization distance threshold 
%                           for pairs of objects. To use the same
%                           threshold among all pairs of objects enter the
%                           same number three times.
%                           colocDistThresh(1): is threshold between 
%                                               channel 1 and 2
%                           colocDistThresh(2): threshold between channel 1 
%                                               and 3
%                           colocDistThresh(3): threshold between channel 2 
%                                               and 3 
%                           Default: [3 3 3] (pixels)
%
%   numTarRefRand:          number of target randomizations to assess conditional
%                           colocalization between target and reference
%                           Default 1:
%
%   numCondRand:            number of condition randomizations to calculate randC
%                           (see vega-Lugo et al. 2022 for randC defintion)
%                           Default: 1
%
%OUPUT Output is saved in directory.   
%      For detailed list of output parameters see condColocGeneral.m documentation
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
%% Input

ip = inputParser;

ip.CaseSensitive = false;
ip.FunctionName  = 'condColocGenralAnalysisMLMD';

%Add required parameter
addRequired(ip, 'MLMD', @(x) isa(x,'MovieData') || isa(x,'MovieList'))
addRequired(ip,'ptOrBlob',@isvector)

%Add optional parameters
addParameter(ip,'detectionProcessID',{'PointSource',[],[]},@iscell)
addParameter(ip,'colocDistThresh', [3 3 3], @isvector)
addParameter(ip,'numTarRefRand',1,@isscalar)
addParameter(ip,'numCondRand',1,@isscalar)

parse(ip,MLMD,ptOrBlob,varargin{:})

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

%go over all movies and run colocalization
for i = 1 : numMovies
    
    %get movieData of current movie
    if listFlag == 1
        MD = MovieData.load(ML.movieDataFile_{i});
    end
    
    iProcCondColoc = MD.getProcessIndex('CondColocGeneralProcess',Inf,0);
    
    if isempty(iProcCondColoc)
       iProcCondColoc = numel(MD.processes_) + 1;
       MD.addProcess(CondColocGeneralProcess(MD));
    end
        
    
    %define colocalization parameters
    
    p = MD.getProcess(iProcCondColoc).getParameters();
    
    p.PtOrBlob = ip.Results.ptOrBlob;
    
    p.DetectionProcessID = ip.Results.detectionProcessID;
        
    p.ColocDistThresh = ip.Results.colocDistThresh;
    p.NumTarRefRand = ip.Results.numTarRefRand;
    p.NumCondRand = ip.Results.numCondRand;
    
    MD.getProcess(iProcCondColoc).setParameters(p);
    MD.getProcess(iProcCondColoc).run
    
end

end