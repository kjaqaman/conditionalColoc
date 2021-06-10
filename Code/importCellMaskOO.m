function importCellMaskOO(processOrMovieData,varargin)
%importCellMask imports a previously-defined cell mask to movieData
%
%Khuloud Jaqaman, March 2015
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

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(processOrMovieData,varargin{:});
paramsIn = ip.Results.paramsIn;

% Use new API to allow the Process to be passed directly
[movieData, maskProc] = getOwnerAndProcess(processOrMovieData,'ImportCellMaskProcess',true);


%Parse input, store in parameter structure
p = parseProcessParams(maskProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check default parameters
if(~isfield(p,'ChannelIndex'))
    p.ChannelIndex = 1;
end
if(~isfield(p,'fileName'))
    p.fileName = cell(size(movieData.channels_));
end
if(~isfield(p,'filePath'))
    p.filePath = cell(size(movieData.channels_));
end
if(~isfield(p,'askUser'))
    p.askUser = true;
end

%Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
mkClrDir(p.OutputDirectory)
for i = p.ChannelIndex;
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
maskProc.setOutFilePaths(outFilePaths);




%% --------------- cell mask import ---------------%%%

for i = p.ChannelIndex
    
    %Ask user where cell mask is currently stored
    if(isempty(p.filePath{i}))
        p.filePath{i} = regexprep(movieData.movieDataPath_,'m(-[0-9][0-9]-01)','BF$1');
    end
    if(isempty(p.fileName{i}))
        p.fileName{i} = 'roiMask.tif';
    end
    fullFileName = [p.filePath{i} filesep p.fileName{i}];
    if(p.askUser || ~exist(fullFileName,'file'))
        [p.fileName{i},p.filePath{i},filterIndx] = uigetfile(fullFileName,['Cell mask file for Channel ' num2str(i) ' of Movie ' movieData.movieDataFileName_(1:end-4)],movieData.roiMaskPath_);
    else
        filterIndx = 1;
    end

    %save cell mask file and its original location
    if filterIndx == 0
        error('No file selected');
    else
        copyfile(fullfile(p.filePath{i},p.fileName{i}),fullfile(p.OutputDirectory,['cellMask_channel_' num2str(i) '.tif']));
        fileName = p.fileName{i};
        filePath = p.filePath{i};
        save(outFilePaths{1,i} ,'fileName','filePath');
    end
    
end

maskProc.setParameters(p);
movieData.save();

