%%---------------------------------------------------------------------
%  SPDX-License-Identifier: BSD-3-Clause
%
%  Copyright(c) 2020 Intel Corporation. All rights reserved.
%
%  Author: Shriram Shastry <malladi.sastry@linux.intel.com>
function drcpowfixed = get_drc_pow(filename)
%---------------------------------------------------
%---------------------------------------
%   History
%---------------------------------------
%   2021/02/09 Sriram Shastry       - initial version
%
%% Initialize variables.

Info = dir(fullfile(filename));
fullFileName            = fullfile(Info.folder,Info.name);
if exist(fullFileName, 'file')
    startRow = 2;
    
    %% Format for each line of text:
    %   column1: double (%f)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    %   column5: double (%f)
    formatSpec = '%11f%12f%12f%12f%f%[^\n\r]';
    
    %% Open the text file.
    fileID = fopen(filename,'r');
    
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    %% Close the text file.
    fclose(fileID);
    
    
    %% Create output variable
    drcpowfixed = table;
    drcpowfixed.idxi = dataArray{:, 1};
    drcpowfixed.idxj = dataArray{:, 2};
    drcpowfixed.inxQ626 = dataArray{:, 3};
    drcpowfixed.inyQ230 = dataArray{:, 4};
    drcpowfixed.outpowQ1220 = dataArray{:, 5};
    
    %% Clear temporary variables
    clearvars filename startRow formatSpec fileID dataArray ans;
else
    % File does not exist.  user warning.
    errorMessage = sprintf('Error: file not found:\n\n%s', fullFileName);
    uiwait(errordlg(errorMessage));
    return;
end
disp(fullFileName);