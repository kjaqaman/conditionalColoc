function plotSaveConditionalColoc(PtOrBlob,plotting, saveRes, varargin)
%PLOTSAVECONDITIONALCOLOC plots and/or saves conditional colocalization results
%
% SYNOPSIS plotSaveConditionalColoc(PtOrBlob,plotting, saveRes, varargin)
%
%INPUT 
%   PtOrBlob:  1 for plotting Pt2Pt2Pt or 2 for Pt2Pt2Blob
%
%   plotting:  1 to plot conditional colocalization results, 0 otherwise
%              Default: 1
%
%   saveRes:   1 to save conditional colocalization results and/or figures, 0 otherwise
%              When saveRes = 1 follow prompt in command window to
%              choose what to save.
%              Default: 0
%           
%   varargin:  MovieList files (minimum 1 file, no maximum)
%              NOTE: function only takes movieList files. 
%
%OUTPUT
% Function outputs plots (when plotting = 1) and compiled colocalization
% measures and number of objects (when saveRes = 1)). When saving files a
% window will show for user to select where to save results.
%
%   Saving colocalizaition results:
%       A mat file with each conditional colocalization measure
%       will be saved in chosen directory. Each each variable
%       (colocalization measure) contains a matrix where:
%
%               rows = number of images
%               Columns contain the following:
%                        Column 1: colocalization measure for original
%                                  positions of reference and target objects
%                        Column 2: colocalization measure after replacing
%                                  target location with grid (nullTR).
%                        Column 3: critical value of colocalization measure
%                                  given input alpha as calculated in
%                                  Helmuth et al. BMC Bioinformatics 2010.
%                                  NOTE: this value is has not been
%                                  validated for rescaled probabilities. It 
%                                  is not used for the analysis in Vega-Lugo
%                                  et al. 2021.
%                        Column 4: ratio of column one to column two
%
%               Third dimesion contains:
%                       First index contains above informaiton for 
%                       original condition object positons
%
%                       Second index contains above information for 
%                       after randomizing condition objects (randC)
%                                             
%
%       Mat file is named colocalMeasureCondARefBTarC where
%       A, B, and C are the index of the channel used for Cond, Ref, and
%       Tar respectively. When A or B say Blob it means non-punctate object.
%       See condColocPt2Pt2Pt and condColocPt2Pt2Blob for .mat file content.
%
%   Saving number of objects:
%       A mat file with number of reference, target, and condition objects as
%       well as cell area will be saved in chosen directory. Number of rows
%       is equal to number of images.
%
%       Mat file is named numObjCondARefBTarC where
%       A, B, and C are the index of the channel used for Cond, Ref, and
%       Tar respectively. When A or B say Blob it means non-punctate object.
%       See condColocPt2Pt2Pt and condColocPt2Pt2Blob for .mat file content.
%
%Jesus Vega-Lugo July 2019
%
%Jesus Vega-Lugo January 2020 generalized the function to plot Pt2Pt2Pt and
%Pt2Pt2Blob. 
%
%Jesus Vega-Lugo December 2020 figures will now show colocalization 
%probability for data and the null probability in the same figure. Made 
%more user friendly by adding an option to save figure, values or both and 
%where to save them
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

%% Parse inputs

%check that type of colocalization to be plotted is specified correctly
if nargin < 1
    error('Especify type of Conditional Coloc. to be plotted. Enter 1 for plotting Pt2Pt2Pt or 2 for PtPt2Blob')
    
elseif PtOrBlob ~= 1 && PtOrBlob ~= 2 
    error('Enter 1 for plotting Pt2Pt2Pt or 2 for PtPt2Blob')
%display text showing the type of colocalization being plotted    
elseif PtOrBlob == 1
    disp('Plotting Pt2Pt2Pt conditional colocalization')
    
elseif PtOrBlob == 2
    disp('Plotting Pt2Pt2Blob conditional colocalization')
end

%check plotting
if nargin < 2 || isempty(plotting)
    plotting = 1;
end

%check save
if nargin < 3 || isempty(saveRes) || saveRes == 0
    whatToSave = 0;
elseif saveRes == 1
    whatToSave = input(['\n' 'Do you want to save the figures, the compiled data values, the number of objects per cell'...
        '\n' 'in a .mat file, or all?'...
        '\n' 'Enter 1 for figures, 2 for compiled data values, 3 number of objects per cell, 4 for all'...
        '\n' '\n' 'Enter number']);
end

%check at least one movie list is input
if isempty(varargin)
    error('At least one MovieList file must be input!!')
end

%% Go through movie lists

%number of movie lists
numMLs = length(varargin);

%number of movies per MovieLists
moviesPerML = NaN(numMLs,1);
%total number movies
totalMovies = 0;


for ml = 1:numMLs
    %check if all varargin inputs are MovieLists
    if isa(varargin{ml},'MovieList')
        moviesPerML(ml) = length(varargin{ml}.movieDataFile_);
        totalMovies = totalMovies + moviesPerML(ml);
        
    else
        %insert error if one of the varargin inputs is not a MovieList
        error('All varargin inputs must be MovieList file')
    end
end

%% Store colocalization values for all movies
colocSummary = cell(totalMovies,6);
numObjectSummary = cell(totalMovies,6);

tempNumMovies = 0;
%go through all movielist
for ml = 1:numMLs
    %go through every movie on the current MovieList 
    for mv = 1:moviesPerML(ml) 
        
        %check if plotting Pt2Pt2Pt or Pt2Pt2Blob
        if PtOrBlob == 1 %Pt2Pt2P
            
            iProc = varargin{ml}.movies_{mv}.getProcessIndex('CondColocPt2Pt2PtProcess',1);
            currentDirC = fullfile(varargin{ml}.movies_{mv}.getProcess(iProc).funParams_.OutputDirectory,...
                'ConditionalColocPt2Pt2Pt','colocalInfoCond.mat');
            
            currentDirN = fullfile(varargin{ml}.movies_{mv}.getProcess(iProc).funParams_.OutputDirectory,...
                'ConditionalColocPt2Pt2Pt','numOfObjects.mat');
            
        elseif PtOrBlob == 2 %Pt2Pt2Blob
            
            iProc = varargin{ml}.movies_{mv}.getProcessIndex('CondColocPt2Pt2BlobProcess',1);
            currentDirC = fullfile(varargin{ml}.movies_{mv}.getProcess(iProc).funParams_.OutputDirectory,...
                'ColocalizationPt2Pt2Blob','colocalInfo.mat');
            
            currentDirN = fullfile(varargin{ml}.movies_{mv}.getProcess(iProc).funParams_.OutputDirectory,...
                'ColocalizationPt2Pt2Blob','numOfObjects.mat');
            
        end
        
        %colocalization values
        tempC = load(currentDirC);
        fieldNamesC = fieldnames(tempC);
        numCategories = length(fieldNamesC);
        
        %number of objects 
        tempN = load(currentDirN);
        fieldNamesN = fieldnames(tempN);
        
        %store colocalization values of each movie
        c = 1;
        for f = 1:numCategories
            if ~isstruct(tempC.(fieldNamesC{f}))
                colocSummary{mv+tempNumMovies,c} = tempC.(fieldNamesC{f});
                
                numObjectSummary{mv+tempNumMovies,c} = tempN.(fieldNamesN{f});
                c = c + 1;
            else
                numCategories = numCategories - 1;
            end
        end 
    end
    tempNumMovies = tempNumMovies + mv;
    
end 

%% Separate categories

pTwRgivenTwC = NaN(totalMovies,4,2);
pTwRgivenTnC = NaN(totalMovies,4,2);
pTwRwCgivenTwC = NaN(totalMovies,4,2);
pTwRnCgivenTnC = NaN(totalMovies,4,2);
pTwRwC = NaN(totalMovies,4,2);
pTwRnC = NaN(totalMovies,4,2);
pTwR = NaN(totalMovies,4,2);
pRwC = NaN(totalMovies,1,2);
pTwC = NaN(totalMovies,1,2);

allTar = NaN(totalMovies,1);
numTwC = NaN(totalMovies,2);
numTnC = NaN(totalMovies,2);
ROIarea = NaN(totalMovies,1);
allRef = NaN(totalMovies,1);
numRwC = NaN(totalMovies,2);
numRnC = NaN(totalMovies,2);
allCond = NaN(totalMovies,1);
refAreas = NaN(totalMovies,1);
condAreas = NaN(totalMovies,1);

%get the scale for y axis
cTcNullVals = vertcat(colocSummary{:});
cTMaxY = max(cTcNullVals(:,1:2),[],'all'); 
 
maxYVals = round(cTMaxY + 0.05,1);
minY = -0.01;

%go through each condition
 for f = 1:numCategories
     %go through each movie
    for mv = 1:totalMovies
        
        pTwRgivenTwC(mv,:,1) = colocSummary{mv,f}(1,:,1);
        pTwRgivenTwC(mv,:,2) = colocSummary{mv,f}(1,:,2);

        pTwRgivenTnC(mv,:,1) = colocSummary{mv,f}(2,:,1);
        pTwRgivenTnC(mv,:,2) = colocSummary{mv,f}(2,:,2);

        pTwRwCgivenTwC(mv,:,1) = colocSummary{mv,f}(3,:,1);
        pTwRwCgivenTwC(mv,:,2) = colocSummary{mv,f}(3,:,2);

        pTwRnCgivenTnC(mv,:,1) = colocSummary{mv,f}(4,:,1);
        pTwRnCgivenTnC(mv,:,2) = colocSummary{mv,f}(4,:,2);

        pTwRwC(mv,:,1) = colocSummary{mv,f}(5,:,1);
        pTwRwC(mv,:,2) = colocSummary{mv,f}(5,:,2);

        pTwRnC(mv,:,1) = colocSummary{mv,f}(6,:,1);
        pTwRnC(mv,:,2) = colocSummary{mv,f}(6,:,2);

        pTwR(mv,:,1) = colocSummary{mv,f}(7,:,1);
        pTwR(mv,:,2) = colocSummary{mv,f}(7,:,2);
        
        pRwC(mv,1,1) = colocSummary{mv,f}(8,1,1);
        pRwC(mv,1,2) = colocSummary{mv,f}(8,1,2);
        
        pTwC(mv,1,1) = colocSummary{mv,f}(9,1,1);
        pTwC(mv,1,2) = colocSummary{mv,f}(9,1,2);
        
        allTar(mv,1) = numObjectSummary{mv,f}.AllTar(1);
        
        numTwC(mv,1) = numObjectSummary{mv,f}.TwC(1);
        numTwC(mv,2) = numObjectSummary{mv,f}.TwC(2);
        
        numTnC(mv,1) = numObjectSummary{mv,f}.TnC(1);
        numTnC(mv,2) = numObjectSummary{mv,f}.TnC(2);
        
        ROIarea(mv,1) = numObjectSummary{mv,f}.ROIarea;
        
        if PtOrBlob == 1
            allRef(mv,1) = numObjectSummary{mv,f}.AllRef(1); 
            
            numRwC(mv,1) = numObjectSummary{mv,f}.RwC(1);
            numRwC(mv,2) = numObjectSummary{mv,f}.RwC(2);
            
            numRnC(mv,1) = numObjectSummary{mv,f}.RnC(1);
            numRwC(mv,2) = numObjectSummary{mv,f}.RnC(2);
            
            allCond(mv,1) = numObjectSummary{mv,f}.AllCond(1);
            
        elseif PtOrBlob == 2  && isfield(numObjectSummary{mv,f},'RefAreas')
            allCond(mv,1) = numObjectSummary{mv,f}.AllCond;
            
            allRef(mv,1) = numObjectSummary{mv,f}.AllRef(1);
            refAreas(mv,1) = mean(numObjectSummary{mv,f}.RefAreas);
            
            
        elseif PtOrBlob == 2 && isfield(numObjectSummary{mv,f},'CondAreas')
            allRef(mv,1) = numObjectSummary{mv,f}.AllRef(1);
            
            allCond(mv,1) = numObjectSummary{mv,f}.AllCond(1);
            condAreas(mv,1) = mean(numObjectSummary{mv,f}.CondAreas);
        end       
    end
    
    %matrix with collocalization probabilities
    pNullBx = [pTwR(:,1:2,1),pTwR(:,1:2,2),...
               pTwRgivenTwC(:,1:2,1),pTwRgivenTwC(:,1:2,2),...
               pTwRgivenTnC(:,1:2,1),pTwRgivenTnC(:,1:2,2),...
               pTwRwC(:,1:2,1),pTwRwC(:,1:2,2),...
               pTwRnC(:,1:2,1),pTwRnC(:,1:2,2),...
               pTwRwCgivenTwC(:,1:2,1),pTwRwCgivenTwC(:,1:2,2),...
               pTwRnCgivenTnC(:,1:2,1),pTwRnCgivenTnC(:,1:2,2)];
  
       
    %matrix with pRC and pTC
    pTCpRCBx = [pTwC(:,:,1),pTwC(:,:,2),pRwC(:,:,1),pRwC(:,:,2)];
    
       bxs = {pNullBx, pTCpRCBx};
        
    %labels for cT and cT/cNull
    labelX = {'p(TwR)','p(TwR|TwC)','p(TwR|TnC)','p(Tw(RwC))','p(Tw(RnC))','p(Tw(RwC)|TwC)','p(Tw(RnC)|TnC)'};
    pLabelX = {'p(TwC)' 'p(RwC)'};
    
%% Plot data results for current condition (f)
    if plotting
       for bx=1:2
           if bx == 1
               %group boxes in cuadruples
               group = repelem((1:7)',totalMovies*4);
               
               %space holder
               sh = repmat({'x1','x2','x3','x4'},totalMovies,7);

               figure, boxplot(bxs{bx}(:),{group, sh(:)},'notch','on','factorgap',...
                   6,'color',[0 0 1;0 0 0;0 0.5 0;0.85 0.325 0.098],'Width',1,'OutlierSize',10,'Symbol','mo')
                h = gca;
                h.XTick = (2.5:5.6:39);
                h.XTickLabel = labelX;
                h.FontWeight = 'bold';
                
                annotation('textbox',[0.2 0.5 0.3 0.3],'string',{'Blue = Data',...
                 'Black = Ref Null','Green = Cond Null','Red = Ref and Cond Null'},'fitboxtotext','on')
             
               %get x coordinates to plot data points on top of boxes
               g = findobj(h,'Tag','Box');
               for i = 1:length(g)
                   center(i) = mean(g(i).XData([1,6]));
               end

               %repaet x value as many times as there are y values on each box
               xScatter = flip(repelem(center',totalMovies));
               center = [];

               hold on

                 blue = [1:totalMovies,totalMovies*4+1:totalMovies*5,...
                   totalMovies*8+1:totalMovies*9,totalMovies*12+1:totalMovies*13,...
                   totalMovies*16+1:totalMovies*17,totalMovies*20+1:totalMovies*21,...
                   totalMovies*24+1:totalMovies*25];

               scatter(xScatter(blue),bxs{bx}(blue)','MarkerEdgeColor',[0 0 1],'LineWidth',0.75)

               black = [totalMovies+1:totalMovies*2,totalMovies*5+1:totalMovies*6,...
                   totalMovies*9+1:totalMovies*10,totalMovies*13+1:totalMovies*14,...
                   totalMovies*17+1:totalMovies*18,totalMovies*21+1:totalMovies*22,...
                   totalMovies*25+1:totalMovies*26];

               scatter(xScatter(black),bxs{bx}(black)','MarkerEdgeColor',[0 0 0],'LineWidth',0.75)

               green = [totalMovies*2+1:totalMovies*3,totalMovies*6+1:totalMovies*7,...
                   totalMovies*10+1:totalMovies*11,totalMovies*14+1:totalMovies*15,...
                   totalMovies*18+1:totalMovies*19,totalMovies*22+1:totalMovies*23,...
                   totalMovies*26+1:totalMovies*27];

               scatter(xScatter(green),bxs{bx}(green)','MarkerEdgeColor',[0 0.5 0],'LineWidth',0.75)

               red = [totalMovies*3+1:totalMovies*4,totalMovies*7+1:totalMovies*8,...
                   totalMovies*11+1:totalMovies*12,totalMovies*15+1:totalMovies*16,...
                   totalMovies*19+1:totalMovies*20,totalMovies*23+1:totalMovies*24,...
                   totalMovies*27+1:totalMovies*28];
       
       scatter(xScatter(red),bxs{bx}(red)','MarkerEdgeColor',[0.85 0.325 0.098],'LineWidth',0.75)
           else
               %group boxes in pairs
               group = repelem((1:2)',totalMovies*2);
               %space holder
               sh = repmat({'x1','x2'},totalMovies,2);

               %plot
               figure, boxplot(bxs{bx}(:),{group, sh(:)},'notch','on','factorgap',...
                   2,'color',[0 0 1;0 0.5 0],'Width',1,'OutlierSize',10,'Symbol','mo')
               h = gca;
               h.XTick = [1.5 3.5];
               h.XTickLabel = pLabelX;
               h.FontWeight = 'bold';
               
               annotation('textbox',[0.2 0.5 0.3 0.3],'string',{'Blue = Data',...
                            'Green = Cond Null'},'fitboxtotext','on')
                        
               %get x coordinates to plot data points on top of boxes
               g = findobj(h,'Tag','Box');
               for i = 1:length(g)
                   center(i) = mean(g(i).XData([1,6]));
               end

               %repaet x value as many times as there are y values on each box
               xScatter = flip(repelem(center',totalMovies));
               center = [];

               hold on

                 blue = [1:totalMovies,totalMovies*2+1:totalMovies*3];

               scatter(xScatter(blue),bxs{bx}(blue)','MarkerEdgeColor',[0 0 1],'LineWidth',0.75)
               
               green = [totalMovies+1:totalMovies*2,totalMovies*3+1:totalMovies*4];

               scatter(xScatter(green),bxs{bx}(green)','MarkerEdgeColor',[0 0.5 0],'LineWidth',0.75)
           end
       
       ylim([minY maxYVals])
       ylabel('Colocalization Probability','fontweight','bold')
       title([fieldNamesC{f}])
       
       hold off 
       pause(1)
       
%% Save results for current condition (f)
       %save figures
       if whatToSave == 1 || whatToSave == 4
           if bx == 1
                [file,path] = uiputfile('*.fig','Choose where to save figure', ['condProbabilities' fieldNamesC{f} '.fig']); 
           elseif bx == 2
               [file,path] = uiputfile('*.fig','Choose where to save figure', ['pTarRefWithCond' fieldNamesC{f} '.fig']);
           end
           
            try
                savefig([path '/' file])
            catch
                continue
            end
       end
        
       end
    end  
    
    %save raw colocalization values for each category 
    if whatToSave == 2 || whatToSave == 4
        [file, path] = uiputfile('*.mat','Choose where to save compiled colocalization results',[fieldNamesC{f} '.mat']);
        
        try
            save([path '/' file], 'pTwRgivenTwC','pTwRgivenTnC','pTwRwCgivenTwC',...
            'pTwRnCgivenTnC','pTwRwC','pTwRnC','pTwR','pRwC','pTwC')
        catch
            continue
        end
    end
    
    if whatToSave == 3 || whatToSave == 4
        if PtOrBlob == 1
            [file, path] = uiputfile('*.mat','Choose where to save number of objects',[fieldNamesN{f} '.mat']);
            try 
                save([path '/' file], 'allTar','numTwC','numTnC','allCond',...
                    'ROIarea','allRef','numRwC','numRnC')
            catch
                continue
            end
            
        elseif PtOrBlob == 2 && isfield(numObjectSummary{mv,f},'RefAreas')
            [file, path] = uiputfile('*.mat','Choose where to save number of objects',[fieldNamesN{f} '.mat']);
            try
                save([path '/' file], 'allTar','numTwC','numTnC','allCond',...
                    'ROIarea','allRef','refAreas')
            catch
                continue
            end
            
        elseif PtOrBlob == 2 && isfield(numObjectSummary{mv,f},'CondAreas')
            [file, path] = uiputfile('*.mat','Choose where to save number of objects',[fieldNamesN{f} '.mat']);
                
            try 
                save([path '/' file], 'allTar','numTwC','numTnC','allCond',...
                    'ROIarea','condAreas','allRef','numRwC','numRnC')
            catch
                continue
            end
        end
    end
    
 end 
end