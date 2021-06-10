classdef ImportCellMaskProcess < DataProcessingProcess
    %Class definition for cell mask import
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
        
    methods(Access = public)   
        
        function obj = ImportCellMaskProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = ImportCellMaskProcess.getName;
                super_args{3} = @importCellMaskOO;                               
                if isempty(funParams)                                       
                    funParams = ImportCellMaskProcess.getDefaultParams(owner,outputDir);                               
                end
                super_args{4} = funParams;                    
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end        
        
    end
    
    methods(Static)
        
        function name =getName()
            name = 'Cell Mask Import';
        end

        function funParams = getDefaultParams(owner,varargin)
            
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:});
            outputDir=ip.Results.outputDir;
            
            % Define default process parameters
            funParams.OutputDirectory = [outputDir  filesep 'ImportedCellMask'];
            funParams.fileName = cell(size(owner.channels_));
            funParams.filePath = cell(size(owner.channels_));
            funParams.ChannelIndex = 1;
            funParams.askUser = true;
        end
        
    end
    
end
