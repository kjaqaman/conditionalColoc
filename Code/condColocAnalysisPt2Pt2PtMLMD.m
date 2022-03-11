function condColocAnalysisPt2Pt2PtMLMD(MLMD,varargin)
%CONDCOLOCANALYSISPT2PT2PT run conditional colocalization analysis for 3 punctate objects. 
%
%SYNOPSIS function condColocAnalysisPt2Pt2PtMLMD(MLMD,varargin)
%
% Function will run conditional colocalization analysis for three punctate
% objects and will make all permutations of target, reference and condition.
% 
%INPUT 
%   Required
%       MLMD:               MovieList or MovieData object for movie(s) to be analyzed
%
%   Optional input arguments as name-value pairs
%       detectionProcessID: Cell array containing desired detection process 
%                           to be used for each channel with punctate objects.
%                           E.g.{DetectProcessForCh1 DetectProcessForCh2 DetectProcessForCh3}
%                           Options: 'SubRes'and 'PointSource'
%                           NOTE: SubRes refers to "Gaussian mixture-model
%                           fitting" in u-Track.
%                           Default: 'SubRes'(for all channels)
%                           i.e. {'SubRes' 'SubRes' 'SubRes'}
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
%   numRandomizations:      number of randomizations to be used for calculating
%                           randC value. (see vega-Lugo et al. 2022 for randC defintion)
%                           Default: 100
%
%       alphaValue:         Alpha value for comparing data to nullTR
%                           (see colocalMeasurePt2Pt)
%                           Default: 0.05
%                           NOTE: for conditional colocalization analysis
%                           as described in Vega-Lugo et al. 2022 this
%                           parameter is not relevant.
%
%OUPUT Output is saved in directory.   
%      For detailed list of output parameters see condColocPt2Pt2Pt.m documentation
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
ip.FunctionName  = 'condColocAnalysisPt2Pt2PtMLMD';

%Add required parameter
addRequired(ip, 'MLMD', @(x) isa(x,'MovieData') || isa(x,'MovieList'))

%Add optional parameters
addParameter(ip,'detectionProcessID',{'SubRes','SubRes','SubRes'}, @iscell)
addParameter(ip,'colocDistThresh', [3 3 3], @isvector)
addParameter(ip,'numRandomizations',100,@isscalar)
addParameter(ip,'alphaValue', 0.05, @isscalar)

parse(ip,MLMD,varargin{:})

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
    
    iProcCondColoc = MD.getProcessIndex('CondColocPt2Pt2PtProcess',Inf,0);
    
    if isempty(iProcCondColoc)
       iProcCondColoc = numel(MD.processes_) + 1;
       MD.addProcess(CondColocPt2Pt2PtProcess(MD));
    end
        
    
    %define colocalization parameters
    
    p = MD.getProcess(iProcCondColoc).getParameters();
    
    p.DetectionProcessID = ip.Results.detectionProcessID;
        
    p.ColocDistThresh = ip.Results.colocDistThresh;
    p.NumRandomizations = ip.Results.numRandomizations;
    p.AlphaValue = ip.Results.alphaValue;
    
    MD.getProcess(iProcCondColoc).setParameters(p);
    MD.getProcess(iProcCondColoc).run
    
end
end
