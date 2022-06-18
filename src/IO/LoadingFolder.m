function [patternNames,DelayTime,numFolder,Np] = LoadingFolder(initial_path,start_num,finish_num)
% loading the folder containing sub-folders for muliti scan(usually named as number)
% copyed & modified from Lingrong and qifengfeng 2022/5/12
% initial_path: initial path of folder
% start_num: the start scan number; positive integer
% finish_num: the finish scan number; positive integer
%
% patternNames: the detail path for every data-pattern; cell
% DelayTime: delay of all patterns; double; unit (10^-4 mm) [numFolder*numPattern]
% numFolder: number of sub-folders that need to be handled; int
% Np: number of patterns in each sub-folder; int

path = uigetdir(initial_path); % path of folder
listing = dir(path); % sub-folders
dirFlags = [listing.isdir];
% Extract only those that are directories.
listing = listing(dirFlags);
filefolders = {listing(3:end).name}; % sub-folders we need
Nsacn = length(filefolders); % number of folders
folderID = zeros(1,Nsacn);
for ii = 1:Nsacn
folderID(ii) = str2double(filefolders{ii});
end
[~,sortID] = sort(folderID); % Sequential index
filefolders = filefolders(sortID);

if nargin>2
    % not valid input
    if start_num>finish_num || finish_num>Nsacn
        error('Not valid input')
    end
    start=start_num;
    finish=finish_num;
end

if nargin==2
    % not valid input
    if start_num>Nsacn
        error('Not valid input')
    end
    start=start_num;
    finish=Nsacn;
end

if nargin==1
    start=1;
    finish=Nsacn;
end

% detect the pattern number in each file folder
% usually to avoid the last incomplete folder in experiment
numFolder=finish-start+1; % the number of folders that need to be handled
Np = zeros(1,numFolder);
DeleteFlag = zeros(1,numFolder);
for ii = 1:numFolder
    subpath = [path,'\',filefolders{ii+start-1}];
    sublisting = dir([subpath,'\*.tiff']);
    patternNames = {sublisting.name};
    Np(ii) = length(patternNames);
    if Np(ii) ~= median(Np(1:ii))
        DeleteFlag(ii) = 1;
    end
end
Np = median(Np); % the number of images in each sub-folder
filefolders=filefolders(start:finish);
filefolders = filefolders(DeleteFlag == 0);
numFolder = length(filefolders);

disp('===all the scans completed===')
disp(['===',num2str(numFolder),' scans with ',num2str(Np),' diffraction patterns in each scan','==='])
disp(['===subfolders from ',filefolders{1},' to ',filefolders{end},'===']);

patternNames=cell(numFolder,Np); % path of all patterns
DelayTime=zeros(numFolder,Np); % delay of all patterns
for ii=1:numFolder
    subpath = [path,'\',filefolders{ii}];
    sublisting = dir([subpath,'\','*.tiff']);
    patternNames(ii,:) = strcat({sublisting.folder},'\',{sublisting.name}); 
    DelayPos={sublisting.name};
    for jj=1:Np
        DelayTime(ii,jj)= str2double(DelayPos{jj}(1:end-4));
    end
end

disp('===pattern name loaded, delay time loaded===')

end