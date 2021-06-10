function statTestResult = condColocStatTest(saveRes)
%STATTESTRESULT takes compiled results from plotSaveConditionalColoc.m and do ranksum test on the results
%
%SYNOPSIS statTestResult = condColocStatTest(saveRes)
% Function does a ranksum test comparing each conditional colocalization  
% measure to its corresponding nulTR and randC and to pTwR. When prompted user
% should select the mat file with the compiled results saved from
% plotSaveConditionalColoc
%
% Function will open a dialog box for user to select the mat file of
% interest. Then, it will show another dialog box asking where to save the output 
%
%INPUT
%   saveRes: 1 for saving results. 0 otherwise
%
%OUTPUT
%   statTestResult: Table containing p-values from ranksum tests
%
%Jesus Vega-Lugo April 2021
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

%% Load mat file
[file,path] = uigetfile('*.mat','Select compiled condColoc results to be tested');

if path == 0
    error('No file was input! Rerun function to input file')
else
    load(fullfile(path,file))
end

%% Do tests

pTwRTest = ranksum(pTwR(:,1,1),pTwR(:,2,1));

try
    pTwRgivenTwCVSpTwR = ranksum(pTwRgivenTwC(:,1,1),pTwR(:,1,1));
    pTwRgivenTwCVSNullpTwRgivenTwC = ranksum(pTwRgivenTwC(:,1,1),pTwRgivenTwC(:,2,1));
    pTwRgivenTwCVSRandpTwRgivenTwC = ranksum(pTwRgivenTwC(:,1,1),pTwRgivenTwC(:,1,2));
catch
    pTwRgivenTwCVSpTwR = nan;
    pTwRgivenTwCVSNullpTwRgivenTwC = nan;
    pTwRgivenTwCVSRandpTwRgivenTwC = nan;
end

try
    pTwRgivenTnCVSpTwR = ranksum(pTwRgivenTnC(:,1,1),pTwR(:,1,1));
    pTwRgivenTnCVSNullpTwRgivenTnC = ranksum(pTwRgivenTnC(:,1,1),pTwRgivenTnC(:,2,1));
    pTwRgivenTnCVSRandpTwRgivenTnC = ranksum(pTwRgivenTnC(:,1,1),pTwRgivenTnC(:,1,2));
catch
    pTwRgivenTnCVSpTwR = nan;
    pTwRgivenTnCVSNullpTwRgivenTnC = nan;
    pTwRgivenTnCVSRandpTwRgivenTnC = nan;
end

try
    pTwRwCVSpTwR = ranksum(pTwRwC(:,1,1),pTwR(:,1,1));
    pTwRwCVSNullpTwRwC = ranksum(pTwRwC(:,1,1),pTwRwC(:,2,1));
    pTwRwCVSRandpTwRwC = ranksum(pTwRwC(:,1,1),pTwRwC(:,1,2));
catch
    pTwRwCVSpTwR = nan;
    pTwRwCVSNullpTwRwC = nan;
    pTwRwCVSRandpTwRwC = nan;
end

try
    pTwRnCVSpTwR = ranksum(pTwRnC(:,1,1),pTwR(:,1,1));
    pTwRnCVSNullpTwRnC = ranksum(pTwRnC(:,1,1),pTwRnC(:,2,1));
    pTwRnCVSRandpTwRnC = ranksum(pTwRnC(:,1,1),pTwRnC(:,1,2));
catch
    pTwRnCVSpTwR = nan;
    pTwRnCVSNullpTwRnC = nan;
    pTwRnCVSRandpTwRnC = nan;
    
end

try
    pTwRwCgivenTwCVSpTwR = ranksum(pTwRwCgivenTwC(:,1,1),pTwR(:,1,1));
    pTwRwCgivenTwCVSNullpTwRwCgivenTwC = ranksum(pTwRwCgivenTwC(:,1,1),pTwRwCgivenTwC(:,2,1));
    pTwRwCgivenTwCSRandpTwRwCgivenTwC = ranksum(pTwRwCgivenTwC(:,1,1),pTwRwCgivenTwC(:,1,2));
catch
    pTwRwCgivenTwCVSpTwR = nan;
    pTwRwCgivenTwCVSNullpTwRwCgivenTwC = nan;
    pTwRwCgivenTwCSRandpTwRwCgivenTwC = nan;
end

try
    pTwRnCgivenTnCVSpTwR = ranksum(pTwRnCgivenTnC(:,1,1),pTwR(:,1,1));
    pTwRnCgivenTnCVSNullpTwRnCgivenTnC = ranksum(pTwRnCgivenTnC(:,1,1),pTwRnCgivenTnC(:,2,1));
    pTwRnCgivenTnCVSRandpTwRnCgivenTnC = ranksum(pTwRnCgivenTnC(:,1,1),pTwRnCgivenTnC(:,1,2));
catch
    pTwRnCgivenTnCVSpTwR = nan;
    pTwRnCgivenTnCVSNullpTwRnCgivenTnC = nan;
    pTwRnCgivenTnCVSRandpTwRnCgivenTnC = nan;
end

try
    pTwCTest = ranksum(pTwC(:,:,1),pTwC(:,:,2));
catch
    pTwCTest = nan;
end

try
    pRwCTest = ranksum(pRwC(:,:,1),pRwC(:,:,2));
catch
    pRwCTest = nan;
end

VSAll = [nan;pTwRgivenTwCVSpTwR;pTwRgivenTnCVSpTwR;pTwRwCVSpTwR;pTwRnCVSpTwR;...
         pTwRwCgivenTwCVSpTwR;pTwRnCgivenTnCVSpTwR;nan;nan];
  
VSNull = [pTwRTest;pTwRgivenTwCVSNullpTwRgivenTwC;pTwRgivenTnCVSNullpTwRgivenTnC;pTwRwCVSNullpTwRwC;...
          pTwRnCVSNullpTwRnC;pTwRwCgivenTwCVSNullpTwRwCgivenTwC;pTwRnCgivenTnCVSNullpTwRnCgivenTnC;nan;nan];
      
VSRand = [nan;pTwRgivenTwCVSRandpTwRgivenTwC;pTwRgivenTnCVSRandpTwRgivenTnC;pTwRwCVSRandpTwRwC;...
          pTwRnCVSRandpTwRnC;pTwRwCgivenTwCSRandpTwRwCgivenTwC;pTwRnCgivenTnCVSRandpTwRnCgivenTnC;pTwCTest;pRwCTest];
  
statTestResult = table(VSAll,VSNull,VSRand,'VariableNames',{'VS_p(TwR)','VS_NullTR','VS_RandC'},...
    'RowNames',{'p(TwR)','p(TwR|TwC)','p(TwR|TnC)','p^rs(Tw(RwC))','p^rs(Tw(RnC))','p^rs(Tw(RwC)|TwC)',...
               'p^rs(Tw(RnC)|TnC)','p(TwC)','p(RwC)'});


%% Save output
if saveRes
    [fileSaved, path] = uiputfile('*.mat','Choose where to save table with test results',['statTest' file(15:end)]);
    
        save([path fileSaved],'statTestResult')
end

end