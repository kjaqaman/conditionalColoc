function plotCondColocPvalues(savePlot,valueRange,fileToLoad,saveResPath)
%PLOTCONDCOLOCPVALUES plots the pvalues from condColocStatTest for each colocalization measure on heatmap
%
%Synopsis: plotCondColocPvalues(savePlot,valueRange)
% Function takes the output of condColocStatTest and plots the biggest
% pvalue among the three comparisons of each colocalization measure (CM).
% Heatmap color will show if median data value is above or below median
% value of its respective comparison for each CM. pTwR is compared to its
% nullTR, pTwC & pRwC are compares to their randC, all conditional measures
% are compared to pTwR.Blue means above, red means below. 
%
%INPUT 
%   savePlot:   1 for saving figure, 0 otherwise
%               Default: 0
%
% Optional:
%   valueRange: vector containing range of values used to define colormap's
%               min and max values.
%               Default: [-15 15]
%
%   fileToLoad: full path of the file to be analyzed or name of file
%               (enetered as type char) if the file is on you current path.
%               NOTE: If empty, a dialog box will show asking for the file
%               to be loaded.
%
%   saveResPath: path where stats table will be saved
%                NOTE: If empty, a dialog box will show asking where to
%                save stats table
%
%OUTPUT
%   The plot. If savePlot = 1, user will be prompted to select the desired
%   path to save plot
%
%Jesus Vega-Lugo September 2021
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

%check savePlot
if nargin < 1 || isempty(savePlot)
    savePlot = 0;
end

%check rangeValues
if nargin < 2 || isempty(valueRange)
    valueRange = [-15 15];
end

if nargin < 3 || isempty(fileToLoad)

    % Load mat file
    [file,path] = uigetfile('*.mat','Select statTestResult to be plotted');

    if path == 0
        error('No file was input! Rerun function to input file')
    else
        load(fullfile(path,file))
    end
    
elseif nargin > 2 && ~isempty(fileToLoad)
    try 
        load(fileToLoad)
        [~,fileName,ext] = fileparts(fileToLoad);
        file = [fileName ext];
    catch
        load(fullfile(cd,fileToLoad))
        file = fileToLoad;
    end
end

%% Extract needed info
pValues = NaN(7,1);

%go over all colocalization measures
for colocMeasure = 1:7 
    
    if colocMeasure == 1%compare pTwR to its nullTR
        abvBel = statTestResult{colocMeasure,4};
        
    elseif any(colocMeasure == (2:5))%comapre conditional measures to pTwR
    abvBel = statTestResult{colocMeasure,2};
        
    elseif any(colocMeasure == [6 7])%compare pTwC and pRwC to its randC
        abvBel = statTestResult{colocMeasure,6};
    end
    
    if abvBel == 1%get a positive value for those above and negative for those below its comaprison
        pValues(colocMeasure) = -log10(max(statTestResult{colocMeasure,1:2:end},[],'omitnan'));
    else
        pValues(colocMeasure) = log10(max(statTestResult{colocMeasure,1:2:end},[],'omitnan'));
    end
    
end

%% Plot
yNames = {'p(TwR)','p(TwR|TwC)','p(TwR|TnC)','prs(Tw(RwC))','prs(Tw(RnC))','p(TwC)','p(RwC)'};
xNames = '-log10(pvalue)';
    
figure,imagesc(pValues);
myColorMap = [[ones(60,1), (0:0.01667:0.99)', (0:0.01667:0.99)']; [(0.99:-0.01667:0)', (0.99:-0.01667:0)', ones(60,1)]];
colorbar
caxis(valueRange)
colormap(myColorMap)

title(file(9:end-4))

ax = gca;

ax.YTickLabel = yNames;
ax.XTickLabel = {'', xNames, ''};
set(ax,'FontWeight','bold')

pValues = abs(pValues);
text(ones(7,1),(1:7)',num2str(pValues,'%.3f'),'HorizontalAlignment','center','FontWeight','bold')

for y = 1.5:1:6.5
    hold on, plot([0.5 1.5],[y y],'k','LineWidth',1)
end

% h = heatmap(xNames,yNames,pValues);
% 
% h.Colormap = [[ones(60,1), (0:0.01667:0.99)', (0:0.01667:0.99)']; [(0.99:-0.01667:0)', (0.99:-0.01667:0)', ones(60,1)]];
%[[ones(6,1) (0:0.1667:1)' (0:0.1667:1)']; [(1:-0.1667:0)' (1:-0.1667:0)' ones(6,1)]; [0 0 1]];
% %maxPvalue = max(pValues,[],'omitnan');
% h.ColorLimits = valueRange;

%% Save figure
if savePlot
    if nargin < 4 || isempty(saveResPath)
        [fileSaved, path] = uiputfile('*.mat','Choose where to save figure with p-values',...
                        ['condColocPvalues' file(9:end-3) 'fig']);

            savefig([path fileSaved])
    
    elseif nargin > 3 && ~isempty(saveResPath) 
        savefig([saveResPath '/condColocPvalues' file(9:end-3) 'fig'])
    end

end