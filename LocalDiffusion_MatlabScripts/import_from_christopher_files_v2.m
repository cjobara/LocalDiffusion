function[dt_theo,nb, t, x, y] = import_from_christopher_files_v2(full_file_name)
delimiter = '\t';
startRow = 3;

formatSpec = '%f%s%s%s%f%s%s%s%f%f%f%f%f%f%[^\n\r]';


fileID = fopen(full_file_name,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
nTracks = dataArray{:, 1};
spaceUnits = dataArray{:, 2};
frameInterval = dataArray{:, 1};
timeUnits = dataArray{:, 4};
generationDateTime = dataArray{:, 5};
from    = dataArray{:, 6};

nSpots  = dataArray{:, 10};
t       = dataArray{:, 11};
x       = dataArray{:, 12};
y       = dataArray{:, 13};
z       = dataArray{:, 14};





dt_theo = frameInterval;
nb      = nSpots;
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;