function [data,dates,series,tcode] = fred_qd2(csv_in,target_names,transformation,stationY)
% =========================================================================
% BREAKDOWN OF THE SCRIPT
% Part 1: Load and label FRED-MD data.
% Part 2: Process data -- transform each series to be stationary and remove
%         outliers.
% -------------------------------------------------------------------------
% AUXILIARY FUNCTIONS
% List of auxiliary functions to be saved in same folder as this script.
%
%   prepare_missing() - transforms series based on given transformation
%       numbers
% -------------------------------------------------------------------------
% NOTES
% Authors: Michael W. McCracken and Serena Ng
% Modify by: Nicolas Hardy
% -------------------------------------------------------------------------
% PARAMETERS TO BE CHANGED
% File name of desired FRED-MD vintage

% Type of transformation performed on each series before factors are
% estimated
%   0 --> no transformation
%   1 --> demean only
%   2 --> demean and standardize
%   3 --> recursively demean and then standardize
DEMEAN=2;

% =========================================================================
% PART 1: LOAD AND LABEL DATA
% Load data from CSV file
dum          = importdata(csv_in,',');
fcode        = dum.data(1,:);
tcode        = dum.data(2,:);
rawdata      = dum.data(3:end,:);
series_names = dum.textdata(1,2:end);  % Cell array of strings

% Define target variable names
%target_names = {'GDPC1', 'GDPCTPI', 'FEDFUNDS'};
% Find positions of the target names in series_names
[is_target, pos_in_target] = ismember(series_names, target_names);

% Get the indices in series_names where target_names are found
idx_target = find(is_target);

% Sort idx_target according to the desired target_names order
[~, desired_order] = ismember(target_names, series_names);
desired_order = desired_order(desired_order > 0);  % Remove not-found entries
idx_first = desired_order;

% Get indices of columns where fcode == 1
 idx_fcode = find(fcode == 1|fcode == 0);

% Combine both sets of indices, ensuring uniqueness
idx_combined = unique([idx_first, idx_fcode], 'stable');
% Reorder: first the target series (in specified order), then the rest
idx_rest = setdiff(idx_combined, idx_first, 'stable');
final_idx = [idx_first, idx_rest];

% Subset and reorder variables
fcode   = fcode(final_idx);
tcode   = tcode(final_idx);
rawdata = rawdata(:, final_idx);
series  = series_names(final_idx);

% Month/year of final observation
final_datevec=datevec(dum.textdata(end,1));
final_month=final_datevec(2);
final_year=final_datevec(1);

% Dates (monthly) are of the form YEAR+MONTH/12
% e.g. March 1970 is represented as 1970+3/12
% Dates go from 1959:01 to final_year:final_month (see above)
dates = (1959+3/12:3/12:final_year+final_month/12)';
% T = number of months in sample
T=size(dates,1);
rawdata=rawdata(1:T,:);

% =========================================================================
% PART 2: PROCESS DATA

% Transform raw data to be stationary using auxiliary function prepare_missing()
% Transformation = 1 avoids over differenciating. Transformation = 0 is the
% same as McCracken

if transformation == 1
    % Replace all 6s with 5s
     tcode(tcode == 6) = 5;
    % Replace all 2s with 1s
     tcode(tcode == 2) = 1;
else
end

% change tcodes for the VAR variables
tcode(1,1:3) = stationY;
yt=prepare_missing(rawdata,tcode);

% Reduce sample to usable dates: remove first two months because some
% series have been first differenced
data=yt(3:T,:);
dates=dates(3:T,:);

%% keep just columns with no NaNs
% cols_to_keep = all(~isnan(data), 1);
% % Keep only those columns
% data   = data(:, cols_to_keep);
% series = series(1,cols_to_keep);
% tcode   = tcode(1,cols_to_keep);


