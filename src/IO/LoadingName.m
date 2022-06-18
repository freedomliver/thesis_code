function [patternNames,DelayTime] = LoadingName(path,filefolders,numFolder,numPattern)
% extracting all pattern name and delay time for preprocessing
% 2020/10/12 qifengfeng
%
% path: the main folder path
% filefolders: sub-folders need to be handled; cell
% numFolder: number of sub-folders that need to be handled; int
% numPattern: number of patterns in each sub-folder; int
%
% patternNames: path of all patterns; cell numFolder*numPattern
% DelayTime: delay of all patterns; double numFolder*numPattern; unit 10^-4 mm

patternNames=cell(numFolder,numPattern); % path of all patterns
DelayTime=zeros(numFolder,numPattern); % delay of all patterns
for ii=1:numFolder
    subpath = [path,'\',filefolders{ii}];
    sublisting = dir([subpath,'\','*.tiff']);
    patternNames(ii,:) = strcat({sublisting.folder},'\',{sublisting.name}); 
    DelayPos={sublisting.name};
    for jj=1:numPattern
        DelayTime(ii,jj)= str2double(DelayPos{jj}(1:end-4));
    end
end

disp('===pattern name loaded, delay time loaded===')

end

