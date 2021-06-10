function simulatedColoc = runColocalizationAnalysisSimulation(simulationType,numTarDetect,varargin)
%RUNCOLOCALIZATIONANALYSISSIMULATION generates simulated data and runs colocalization analysis 
%
%SYNOPSIS simulatedColoc = runColocalizationAnalysisSimulation(simulationType,numTarDetect,varargin)
%
% Function compiles input parameters needed to set up simulations and
% subsequent analysis. It outputs colocalization results from simulation. 
% When multiple values are input for the same parameter function will loop to
% each value making all possible parameter combinations given the input parameters.
%
%INPUT
% Needed parameters for all types of simulation
%   simulationType: type of colocalization you want to simulate
%                   Options: 'Pt2Pt', 'Pt2Blob', 'Pt2Pt2Pt', 'Pt2Pt2Blob'
%
%   numTarDetect: scalar or vector containing number of traget objects
%
%%%%Every parameter below this line must be input as a name-value pair%%%%
%
%   colocDistRange: range of distances (in pixels) target will be
%                   positioned away from reference or condition objects
%                   Default: [0 distThresh-2] when simulationType = Pt2Pt 
%                            or Pt2Blob
%                            [0 max(colocDistThresh)-2] when simulationType =
%                            Pt2Pt2 or Pt2Pt2Blob
%                   NOTE:simulation code add noise to the target positions
%                   using rand so, max(colocDistRange) will be ~1 higher
%                   than what is input. E.g. if colocDistRange = [0 1] 
%                   objects will be positioned between 0 and 2.
%
%%%%%%%%%%%%%%Pt2Pt simulation (two punctate objects)%%%%%%%%%%%%%%%%%%
%
% Needed for Pt2Pt simulation
%   numRefDetect: scalar or vector containing number of reference objects
%
%   colocFraction: scalar or vector containing fraction of target objects 
%                  colocalizing with reference.
%   
%   For defining simulation area, either ROIs from data or square's side length
%   can be input:
%       cellMasks: cell array containing binary images of ROIs
%
%       simAreaSideLength: desired length of a square edge.
%                          Example: if simAreaSideLength = 200 the
%                          simulation will be done in a 200 by 200 square
%
%   repeats: desired number of simulations
%            NOTE: this parameter is REQUIRED when using simAreaSideLength.
%            When using cellMasks this parameter is not needed and the number
%            of simulations will be equal to the number of input ROIs.
%
%   distThresh: distance threshold for colocalization analysis
%               Default: 3
%
%%%%%%%%%%%Pt2Blob simulation (1 punctate 1 non-punctate object)%%%%%%%%%%%
%
% Needed parameters for Pt2Blob simulation   
%   segmentations: cell array containing binary images of segmentations of
%                  non-punctate objects. Number of simulations will be equal to
%                  number of segmented images.
%
%       cellMasks: cell array containing binary images of ROIs
%                  NOTE: number of cell masks must be equal to the number
%                  of segmentations. 
%
%   colocFraction: scalar or vector containing fraction of target objects 
%                  colocalizing with non-punctate objects.
%
%   distThresh: distance threshold for colocalization analysis
%               Default: 3
%   
%%%%%%%%%%%%%%Pt2Pt2Pt simulation (three punctate objects)%%%%%%%%%%%%%%%%%
%
% Needed parameters for Pt2Pt2Pt simulation 
%   numRefDetect: scalar or vector containing number of reference objects
%
%   numCondDetect: scalar or vector containing number of conditition objects
%
%   For defining simulation area, either ROIs from data or square's side length
%   can be input:
%       cellMasks: cell array containing binary images of ROIs
%
%       simAreaSideLength: desired length of a square edge.
%                          Example: if simAreaSideLength = 200 the
%                          simulation will be done in a 200 by 200 square
%
%   repeats: desired number of simulations
%            NOTE: this parameter is REQUIRED when using simAreaSideLength.
%            When using cellMasks this parameter is not needed and the number
%            of simulations will be equal to the number of input ROIs.
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
%   colocDistThresh: Vector containing distance thresholds for conditional
%                    colocalization analysis
%                    colocDistThresh(1): threshold between 
%                                        reference and target
%                    colocDistThresh(2): threshold between reference
%                                        and condition
%                    colocDistThresh(3): threshold between target
%                                        and condition
%                    Defualt: [3 3 3] (pixels)
%
%%%%%%%%%%Pt2Pt2Blob simulation (two punctate, one non-punctate)%%%%%%%%%%%
%
% Needed parameters for Pt2Pt2Blob simulation
%   segmentations: cell array containing binary images of segmentations of
%                  non-punctate objects. Number of simulations will be equal to
%                  number of segmented images.
%
%       cellMasks: cell array containing binary images of ROIs
%                  NOTE: number of cell masks must be equal to the number
%                  of segmentations. 
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
%   colocDistThresh: Vector containing distance thresholds for conditional
%                    colocalization analysis
%                    colocDistThresh(1): threshold between 
%                                        reference and target
%                    colocDistThresh(2): threshold between reference
%                                        and condition
%                    colocDistThresh(3): threshold between target
%                                        and condition
%                    Defualt: [3 3 3] (pixels)
%
%OUTPUT
%   simulatedColoc: structure array with number of entries equal to the
%                   number of parameter combinations, containing the
%                   following fields:
%
%                   Fields containing results: 
%                       colocalization: contains a cell array where each
%                                       row represents one simulation. See
%                                       core function of each type of
%                                       simulation for more details on
%                                       colocalization output.
%
%                       numOfObjects: contains cell array where rows
%                                     represents one simulation. See core
%                                     function of each type of simulation
%                                     for more details on colocalization
%                                     output. NOTE: this field is only
%                                     relevant for pt2pt2pt and pt2pt2blob
%                                     simulations.
%
%                   Remaining fields describe input parameters as relevant
%                   for each type of simulation. Field names are identical
%                   to input paramters names.
%
%Jesus Vega-Lugo Dec. 2020
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

%% Parse Parameters

ip = inputParser;
ip.CaseSensitive = false;
ip.FunctionName = 'runColocalizationSimulation';

addRequired(ip,'simulationType',@ischar);
addRequired(ip,'numTarDetect',@(x) isscalar(x) || isvector(x))

addParameter(ip,'colocDistRange',[],@isvector)
addParameter(ip,'numRefDetect',@(x) isscalar(x) || isvector(x))
addParameter(ip,'numCondDetect',[], @(x) isscalar(x) || isvector(x))

addParameter(ip,'segmentations',[],@(x) iscell(x) || islogical(x))
addParameter(ip,'cellMasks',[],@(x) iscell(x) || islogical(x))

addParameter(ip,'simAreaSideLength',[],@isscalar)
addParameter(ip,'repeats',@isscalar)

addParameter(ip,'colocFraction',[],@(x) isscalar(x) || isvector(x))
addParameter(ip,'refColocWithCondFrac',[], @(x) isscalar(x) || isvector(x))
addParameter(ip,'tarColocWithCondFrac',[], @(x) isscalar(x) || isvector(x))
addParameter(ip,'tarColocWithRefAtCondFrac',[], @(x) isscalar(x) || isvector(x))
addParameter(ip,'tarColocWithRefNotAtCondFrac',[], @(x) isscalar(x) || isvector(x))

addParameter(ip,'distThresh',3,@isscalar)
addParameter(ip,'colocDistThresh',[3 3 3],@isvector)

parse(ip,simulationType,numTarDetect,varargin{:})

simulationType = ip.Results.simulationType;
numTarDetect = ip.Results.numTarDetect;
colocDistRange = ip.Results.colocDistRange;

numRefDetect = ip.Results.numRefDetect;
numCondDetect= ip.Results.numCondDetect;

segmentations = ip.Results.segmentations;
cellMasks = ip.Results.cellMasks;

simAreaSideLength = ip.Results.simAreaSideLength;
repeats = ip.Results.repeats;

colocFraction = ip.Results.colocFraction;

refColocWithCondFrac = ip.Results.refColocWithCondFrac;
tarColocWithCondFrac = ip.Results.tarColocWithCondFrac;
tarColocWithRefAtCondFrac = ip.Results.tarColocWithRefAtCondFrac;
tarColocWithRefNotAtCondFrac = ip.Results.tarColocWithRefNotAtCondFrac;

distThresh = ip.Results.distThresh;
colocDistThresh = ip.Results.colocDistThresh;

if isempty(colocDistRange)
    switch simulationType

        case {'pt2pt' 'Pt2Pt' 'Pt2pt' 'pt2Pt' 'pt2blob' 'Pt2Blob' 'Pt2blob' 'pt2Blob'}
            colocDistRange = [0 distThresh - 2];

        case {'pt2pt2pt' 'Pt2Pt2Pt' 'pt2pt2blob' 'Pt2Pt2Blob'}
            colocDistRange =[0  max(colocDistThresh) - 2];
            
        otherwise
            error('Valid simulationType must be input. See function documentaiton.')
    end
end

%% Set up simulations 

switch simulationType
    
    case {'pt2pt' 'Pt2Pt' 'Pt2pt' 'pt2Pt'}
        
        c = 0;
        %get number of values per parameter
        valsNumRefDetect = length(numRefDetect);
        valsNumTarDetect = length(numTarDetect);
        vecColocFrac = length(colocFraction);
        
        if valsNumRefDetect == 0 
            error('Must input at least one value for numRefDetect.')
            
        elseif valsNumTarDetect == 0 
            error('Must input at least one value for numTarDetect.')
            
        elseif vecColocFrac == 0
            error('Must input at least one value for colocFraction.')
        end
        
        %prealocate space for output
        simulatedColoc = repmat(struct('simulationType',[],'numRefDetect',[],'numTarDetect',[],...
            'colocDistRange',[],'colocFraction',[],'repeats',[],'distThresh',[],...
            'colocalization',[]),valsNumRefDetect*valsNumTarDetect*vecColocFrac,1);
        
        %go over ref values
        for iRef = 1:valsNumRefDetect
            
            if ~isempty(cellMasks) %when masks are input use them
                repeats = length(cellMasks); 

                ptMasks = cell(repeats,1);
                for i = 1:repeats
                    img = zeros(size(cellMasks{i}));
                    
                    erodedMask = cellMasks{i};
                    erodedMask(1,:) = 0;
                    erodedMask(size(erodedMask,1),:) = 0;
                    erodedMask(:,1) = 0;
                    erodedMask(:,size(erodedMask,2)) = 0;
                    
                    %erosion to avoid condition detection being to close to
                    %cell edge
                    erodedMask = imerode(erodedMask,strel('square',10));
                    maskIdxs = find(erodedMask);
                    
                    condPositions = datasample(maskIdxs,numRefDetect(iRef),'Replace',false);

                    img(condPositions) = 1;
                    ptMasks{i} = logical(img);
                end
                usingSquare = 0;
            else %when no mask input create area for simulation
                
                if isempty(repeats) %show error if no mask and no repeats or simAreaSideLength input
                    error('You must enter the number of repeats if no cell masks are input')
                elseif isempty(simAreaSideLength)
                    error('You must enter simAreaSideLength when cell masks are not input')
                end
                
                ptMasks = cell(repeats,1);
                cellMasks = cell(repeats,1);
                
                %create simulation area
                maskSz = simAreaSideLength + 100;
                mask = zeros(maskSz);
                mask(maskSz/2,maskSz/2) = 1;
                mask = imdilate(mask,strel('square',simAreaSideLength));
                
                for i = 1:repeats
                    img = zeros(maskSz);
                    
                    %erosion to avoid condition detection being to close to
                    %cell edge
                    erodedMask = imerode(mask,strel('square',10));
                    maskIdxs = find(erodedMask);
                    
                    condPositions = datasample(maskIdxs,numRefDetect(iRef),'Replace',false);

                    img(condPositions) = 1;
                    ptMasks{i} = logical(img);%mask with ref positions
                    cellMasks{i} = logical(mask);
                end
                %falg to set output field simAreaSideLength
                usingSquare = 1;
            end
            
            % go over target values
            for iTar = 1:valsNumTarDetect
                %go over colocalization fraction values
                for iColoc = 1:vecColocFrac
                    
                    colocalization = NaN(repeats,3);
                    
                    %repeat simulation 
                    for iRepeat = 1:repeats
                        [colocalization(iRepeat,1),colocalization(iRepeat,2),...
                            colocalization(iRepeat,3)] = simulatePt2BlobColoc(numTarDetect(iTar),...
                            ptMasks{iRepeat},cellMasks{iRepeat},colocFraction(iColoc),...
                            'colocDistRange',colocDistRange,'distThresh',distThresh);
                    end
                    
                    %store simulation parameters and colocalization results
                    c = c + 1;
                    simulatedColoc(c).simulationType = simulationType;
                    simulatedColoc(c).numRefDetect = numRefDetect(iRef);
                    simulatedColoc(c).numTarDetect = numTarDetect(iTar);
                    simulatedColoc(c).colocDistRange = colocDistRange;
                    simulatedColoc(c).colocFraction = colocFraction(iColoc);
                    simulatedColoc(c).repeats = repeats;
                    simulatedColoc(c).distThresh = distThresh;
                    simulatedColoc(c).colocalization = colocalization;
                    if usingSquare == 1
                        simulatedColoc(c).simAreaSideLength = simAreaSideLength;
                    end
                    
                end
            end
        end
        
    case{'pt2blob' 'Pt2Blob' 'Pt2blob' 'pt2Blob'}
        
        repeats = length(segmentations);
        if repeats == 0
            error('A cell array with segmentations of non-punctate objects must be input. See function documentation.')
        
        elseif isempty(cellMasks) || length(cellMasks) ~= repeats
            error('A cell array with the same number of cell masks as there segmentaitons must be input. See function documentation.')
            
        end
        %get number of values per parameter
        valsNumTarDetect = length(numTarDetect);
        vecColocFrac = length(colocFraction);
        
        if valsNumTarDetect == 0
            error('Must input at least one value for numTarDetect.')
            
        elseif vecColocFrac == 0
            error('Must input at least one value for colocFraction')
        end
        
        %prealocate space for output
        simulatedColoc = repmat(struct('simulationType',[],'numTarDetect',[],'colocDistRange',[],...
            'colocFraction',[],'repeats',[],'distThresh',[],'colocalization',[]),...
            valsNumTarDetect*vecColocFrac,1);
        
        c = 0;
        %go over target values
        for iTar = 1:valsNumTarDetect
            %go over colocalization fraction values
            for iColocFrac = 1:vecColocFrac
                colocalization = NaN(repeats,3);
                
                %repeat simulation
                for iSeg = 1:repeats
                    [colocalization(iSeg,1),colocalization(iSeg,2),colocalization(iSeg,3)] =...
                        simulatePt2BlobColoc(numTarDetect(iTar),segmentations{iSeg},...
                        cellMasks{iSeg},colocFraction(iColocFrac),'colocDistRange',...
                        colocDistRange,'distThresh',3);
                end
                %store simulation parameters and colocalization results
                c = c + 1;
                simulatedColoc(c).simulationType = simulationType;
                simulatedColoc(c).numTarDetect = numTarDetect(iTar);
                simulatedColoc(c).colocDistRange = colocDistRange;
                simulatedColoc(c).colocFraction = colocFraction(iColocFrac);
                simulatedColoc(c).repeats = repeats;
                simulatedColoc(c).distThresh = distThresh;
                simulatedColoc(c).colocalization = colocalization;
            end
        end
        
    case {'pt2pt2pt' 'Pt2Pt2Pt'}
        
        %get number of values per parameter
        valsNumRefDetect = length(numRefDetect);
        valsNumTarDetect = length(numTarDetect);
        valsNumCondDetect = length(numCondDetect);
        
        valsRefColocWithCondFrac = length(refColocWithCondFrac);
        valsTarColocWithCondFrac = length(tarColocWithCondFrac);
        valsTarColocWithRefAtCondFrac = length(tarColocWithRefAtCondFrac);
        valsTarColocWithRefNotAtCondFrac = length(tarColocWithRefNotAtCondFrac);
        
        if valsNumRefDetect == 0 
            error('Must input at least one value for numRefDetect.')
            
        elseif valsNumTarDetect == 0 
            error('Must input at least one value for numTarDetect.')
            
        elseif valsNumCondDetect == 0
            error('Must input at least one value for numCondDetect.')
            
        elseif valsRefColocWithCondFrac == 0
            error('Must input at least one value for valsRefColocWithCondFrac')
            
        elseif valsTarColocWithCondFrac == 0
            error('Must input at least one value for valsTarColocWithCondFrac')
            
        elseif valsTarColocWithRefAtCondFrac == 0
            error('Must input at least one value for valsTarColocWithRefAtCondFrac')
            
        elseif valsTarColocWithRefNotAtCondFrac == 0
            error('Must input at least one value for valsTarColocWithRefNotAtCondFrac')
            
        end
        
        %get number of parameter combinations
        combinations = valsNumRefDetect*valsNumTarDetect*valsNumCondDetect*...
            valsRefColocWithCondFrac*valsTarColocWithCondFrac*valsTarColocWithRefAtCondFrac*...
            valsTarColocWithRefNotAtCondFrac;
        
        simulatedColoc = repmat(struct('simulationType',[],'numRefDetect',[],...
            'numTarDetect',[],'numCondDetect',[],'colocDistRange',[],'refColocWithCond',[],...
            'tarColocWithCond',[],'tarColocWithRefAtCond',[],...
            'tarColocWithRefNotAtCond',[],'repeats',[],'colocDistThresh',[],...
            'colocalization',[],'numOfObjects',[]),combinations,1);
        
        c = 0;
        %go over ref values
        for iRef = 1:valsNumRefDetect
        %go over target values
        for iTar = 1:valsNumTarDetect
        %go over condition values    
        for iCond = 1:valsNumCondDetect
            
            if ~isempty(cellMasks)%when masks are inut use them
                repeats = length(cellMasks);

                ptMasks = cell(repeats,1);
                for i = 1:repeats
                    img = zeros(512);
                    
                    erodedMask = cellMasks{i};
                    erodedMask(1,:) = 0;
                    erodedMask(size(erodedMask,1),:) = 0;
                    erodedMask(:,1) = 0;
                    erodedMask(:,size(erodedMask,2)) = 0;
                    
                    %erosion to avoid condition detection being to close to
                    %cell edge
                    erodedMask = imerode(erodedMask,strel('square',10));
                    maskIdxs = find(erodedMask);
                    
                    condPositions = datasample(maskIdxs,numCondDetect(iCond),'Replace',false);

                    img(condPositions) = 1;
                    ptMasks{i} = logical(img);
                end
                usingSquare = 0;
            else %when no mask input create area for simulation
                
                if isempty(repeats)%show error if no mask and no repeats or simAreaSideLength input
                    error('You must enter the number of repeats if no cell masks are input')
                elseif isempty(simAreaSideLength)
                    error('You must enter simArea when cell masks are not input')
                end
                
                ptMasks = cell(repeats,1);
                cellMasks = cell(repeats,1);
                
                %create simulation area
                maskSz = simAreaSideLength + 100;
                mask = zeros(maskSz);
                mask(maskSz/2,maskSz/2) = 1;
                mask = imdilate(mask,strel('square',simAreaSideLength));
                
                for i = 1:repeats
                    img = zeros(maskSz);
                    
                    %erosion to avoid condition detection being to close to
                    %cell edge
                    erodedMask = imerode(mask,strel('square',10));
                    maskIdxs = find(erodedMask);
                    
                    condPositions = datasample(maskIdxs,numCondDetect(iCond),'Replace',false);

                    img(condPositions) = 1;
                    ptMasks{i} = logical(img);
                    cellMasks{i} = logical(mask);
                end
                %flag to set output field simAreaSideLength
                usingSquare = 1;
            end
                %go over all colocalization fraction values
                for iColocRefWithCond = 1:valsRefColocWithCondFrac
                for iColocTarWithCond = 1:valsTarColocWithCondFrac
                for iColocTarWithRefAtCond = 1:valsTarColocWithRefAtCondFrac
                for iColocTarWithRefNotAtCond = 1:valsTarColocWithRefNotAtCondFrac
                    
                    colocMeasure = cell(repeats,1);
                    numOfObjects = cell(repeats,1);
                    r = 1;
                    %repeat simulation
                    for iRepeat = 1:repeats
                        [colocMeasure{r,1}, numOfObjects{r,1}]= simulatePt2Pt2BlobCondColoc(numRefDetect(iRef),numTarDetect(iTar),...
                            ptMasks{iRepeat},cellMasks{iRepeat},refColocWithCondFrac(iColocRefWithCond),...
                            tarColocWithCondFrac(iColocTarWithCond),tarColocWithRefAtCondFrac(iColocTarWithRefAtCond),...
                            tarColocWithRefNotAtCondFrac(iColocTarWithRefNotAtCond),...
                            'colocDistRange',colocDistRange,'colocDistThresh',colocDistThresh);
                        r = r + 1;
                    end
                    %store simulation parameters and colocalization results
                    c = c +1;
                    simulatedColoc(c).simulationType = simulationType;
                    simulatedColoc(c).numRefDetect = numRefDetect(iRef);
                    simulatedColoc(c).numTarDetect = numTarDetect(iTar);
                    simulatedColoc(c).numCondDetect = numCondDetect(iCond);
                    simulatedColoc(c).colocDistRange = colocDistRange;
                    simulatedColoc(c).refColocWithCond = refColocWithCondFrac(iColocRefWithCond);
                    simulatedColoc(c).tarColocWithCond = tarColocWithCondFrac(iColocTarWithCond);
                    simulatedColoc(c).tarColocWithRefAtCond = tarColocWithRefAtCondFrac(iColocTarWithRefAtCond);
                    simulatedColoc(c).tarColocWithRefNotAtCond = tarColocWithRefNotAtCondFrac(iColocTarWithRefNotAtCond);
                    simulatedColoc(c).repeats = repeats;
                    simulatedColoc(c).colocDistThresh = colocDistThresh;
                    simulatedColoc(c).colocalization = colocMeasure;
                    simulatedColoc(c).numOfObjects = numOfObjects;
                    if usingSquare == 1
                        simulatedColoc(c).simAreaSideLength = simAreaSideLength;
                    end
                    
                end
                end
                end
                end
        end
        end
        end
        
    case {'pt2pt2blob' 'Pt2Pt2Blob'}
        
        %get number of values per parameter
        valsNumRefDetect = length(numRefDetect);
        valsNumTarDetect = length(numTarDetect);
        
        valsRefColocWithCondFrac = length(refColocWithCondFrac);
        valsTarColocWithCondFrac = length(tarColocWithCondFrac);
        valsTarColocWithRefAtCondFrac = length(tarColocWithRefAtCondFrac);
        valsTarColocWithRefNotAtCondFrac = length(tarColocWithRefNotAtCondFrac);
        
        if valsNumRefDetect == 0 
            error('Must input at least one value for numRefDetect.')
            
        elseif valsNumTarDetect == 0 
            error('Must input at least one value for numTarDetect.')
            
        elseif valsRefColocWithCondFrac == 0
            error('Must input at least one value for valsRefColocWithCondFrac')
            
        elseif valsTarColocWithCondFrac == 0
            error('Must input at least one value for valsTarColocWithCondFrac')
            
        elseif valsTarColocWithRefAtCondFrac == 0
            error('Must input at least one value for valsTarColocWithRefAtCondFrac')
            
        elseif valsTarColocWithRefNotAtCondFrac == 0
            error('Must input at least one value for valsTarColocWithRefNotAtCondFrac')
            
        end
        
        %number of parameter combinations
        combinations = valsNumRefDetect*valsNumTarDetect*...
            valsRefColocWithCondFrac*valsTarColocWithCondFrac*valsTarColocWithRefAtCondFrac*...
            valsTarColocWithRefNotAtCondFrac;
        
        %prealocate space for output
        simulatedColoc = repmat(struct('simulationType',[],'numRefDetect',[],'numTarDetect',[],...
            'colocDistRange',[],'refColocWithCond',[],'tarColocWithCond',[],...
            'tarColocWithRefAtCond',[],'tarColocWithRefNotAtCond',[],'repeats',[],...
            'colocDistThresh',[],'colocalization',[],'numOfObjects',[]),combinations,1);
        
        repeats = length(segmentations);
        
        if repeats == 0
            error('A cell array with segmentations of non-punctate objects must be input. See function documentation.')
        
        elseif isempty(cellMasks) || length(cellMasks) ~= repeats
            error('A cell array with the same number of cell masks as there segmentaitons must be input. See function documentation.')
            
        end
        
        c = 0;
        %go over reference values
        for iRef = 1:valsNumRefDetect
        %go over target values
        for iTar = 1:valsNumTarDetect
                %go over colocalization fractions values
                for iColocRefWithCond = 1:valsRefColocWithCondFrac
                for iColocTarWithCond = 1:valsTarColocWithCondFrac
                for iColocTarWithRefAtCond = 1:valsTarColocWithRefAtCondFrac
                for iColocTarWithRefNotAtCond = 1:valsTarColocWithRefNotAtCondFrac
                    
                    colocMeasure = cell(repeats,1);
                    numOfObjects = cell(repeats,1);
                    r = 1;
                    
                    %repeat simulation
                    for iRepeat = 1:repeats
                        [colocMeasure{r,1}, numOfObjects{r,1}] = simulatePt2Pt2BlobCondColoc(numRefDetect(iRef),numTarDetect(iTar),...
                            segmentations{iRepeat},cellMasks{iRepeat},refColocWithCondFrac(iColocRefWithCond),...
                            tarColocWithCondFrac(iColocTarWithCond),tarColocWithRefAtCondFrac(iColocTarWithRefAtCond),...
                            tarColocWithRefNotAtCondFrac(iColocTarWithRefNotAtCond),...
                            'colocDistRange',colocDistRange,'colocDistThresh',colocDistThresh);
                        r = r + 1;
                    end
                    %store simulation parameters and colocalization results
                    c = c +1;
                    simulatedColoc(c).simulationType = simulationType;
                    simulatedColoc(c).numRefDetect = numRefDetect(iRef);
                    simulatedColoc(c).numTarDetect = numTarDetect(iTar);
                    simulatedColoc(c).colocDistRange = colocDistRange;
                    simulatedColoc(c).refColocWithCond = refColocWithCondFrac(iColocRefWithCond);
                    simulatedColoc(c).tarColocWithCond = tarColocWithCondFrac(iColocTarWithCond);
                    simulatedColoc(c).tarColocWithRefAtCond = tarColocWithRefAtCondFrac(iColocTarWithRefAtCond);
                    simulatedColoc(c).tarColocWithRefNotAtCond = tarColocWithRefNotAtCondFrac(iColocTarWithRefNotAtCond);
                    simulatedColoc(c).repeats = repeats;
                    simulatedColoc(c).colocDistThresh = colocDistThresh;
                    simulatedColoc(c).colocalization = colocMeasure;
                    simulatedColoc(c).numOfObjects = numOfObjects;
                end
                end
                end
                end
        end
        end
        
  otherwise
       error('Valid simulationType must be input. See function documentaiton.')
end

end