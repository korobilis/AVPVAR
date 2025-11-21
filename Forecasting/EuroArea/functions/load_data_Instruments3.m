function   [Y,Z,T,M,data,dataZ,dates,varnames] = load_data_Instruments3(select,selectZ,tcode,tcodeZ,standard,standardZ)
%% Function to load data for VAR
% read in VAR series
C          = readtable('data.xlsx');


All_data  = table2array(C(:,2:end));
varnames   = C.Properties.VariableNames(2:end);
dates      = table2array(C(:,1));
[~, index] = ismember(select,varnames); 
data       = All_data(:,index); %VAR data 

[~, indexZ] = ismember(selectZ,varnames); 
dataZ      = All_data(:,indexZ); %Instruments data

% Transform data to stationarity
Y = NaN(size(data));
for i = 1:size(data,2)
    if tcode(i) == 1      % levels
        Y(:,i) = data(:,i);
    elseif tcode(i) == 2  % first differences
        Y(:,i) = [NaN; 100*diff(data(:,i))];
    elseif tcode(i) == 4  % log-levels
        Y(:,i) = log(data(:,i));
    elseif tcode(i) == 5  % first log differences
        Y(:,i) = [NaN; 100*diff(log(data(:,i)))];
    end
end
Y = Y(2:end,:);

if standard == 1
    Y = normalize(Y);
end

% Transform Instruments to stationarity
Z = NaN(size(dataZ));
for i = 1:size(dataZ,2)
    if tcodeZ(i) == 1      % levels
        Z(:,i) = dataZ(:,i);
    elseif tcodeZ(i) == 2  % first differences
        Z(:,i) = [NaN; 100*diff(dataZ(:,i))];
    elseif tcodeZ(i) == 4  % log-levels
        Z(:,i) = log(dataZ(:,i));
    elseif tcodeZ(i) == 5  % first log differences
        Z(:,i) = [NaN; 100*diff(log(dataZ(:,i)))];
    end
end
Z = Z(2:end,:);

if standardZ == 1
    Z = normalize(Z);
end

%% ==== Keep sample periods with no NaNs ===
% Find the first row where all elements are numbers (no NaN)
rowIndex = find(all(~isnan([Y Z]), 2), 1, 'first');
Y = Y(rowIndex:end,:);
Z = Z(rowIndex:end,:);
dates = dates(rowIndex:end);

% Keep just the rows where all variables Y and Z are numbers (discard all
% the rows if some variables end with NaN)
% Find the last row where all elements are numbers
lastValidRow = find(all(~isnan([Y Z]), 2), 1, 'last');
Y = Y(1:lastValidRow, :);
Z = Z(1:lastValidRow, :);
dates = dates(1:lastValidRow);

[T,M] = size(Y);
% [~,Mz] = size(Z);



