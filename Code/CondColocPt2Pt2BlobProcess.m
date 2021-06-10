classdef CondColocPt2Pt2BlobProcess < ImageAnalysisProcess
   % A concreate class for measuring conditional colocalization between 3 images
   % Anthony Vega 09/2014
   %Adapted by Jesus Vega-Lugo (May, 2019) to run conditional colocalization 
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
        methods (Access = public)
        function obj = CondColocPt2Pt2BlobProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData') || isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = CondColocPt2Pt2BlobProcess.getName;
                super_args{3} = @condColocPt2Pt2BlobWrapper;
                if isempty(funParams)
                    funParams = CondColocPt2Pt2BlobProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                if(nargin > 4)
                    super_args{5:nargin} = varargin{5:nargin};
                end
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'CondColocPt2Pt2Blob';
        end

        function methods = getMethods(varargin)
            colocalizationMethods(1).name = 'Point2Point2Blob';
            colocalizationMethods(1).func = @condColocPt2Pt2Blob;            
            
            ip=inputParser;
            ip.addOptional('index',1:length(colocalizationMethods),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            methods=colocalizationMethods(index);
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData') || isa(x,'MovieList'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.DetectedChannels = [1 3];
            funParams.ChannelBlob = 2;
            funParams.ColocDistThresh = [3 3 3];
            funParams.AlphaValue = 0.05;
            funParams.OutputDirectory = [outputDir  filesep];
            funParams.ProcessIndex = [];
            funParams.DetectionProcessID = {'SubRes', 'SubRes'};
        end
    end
end