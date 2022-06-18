function [outMaskpattern,numRegionRemoved] = HotPixelRemoval(inPattern,inMaskPattern,RegionRemoved)
% removing hot piexl caused by X-Ray at EMCCD
% 2020/10/15 qifengfeng
%
% inPattern: input pattern; 2D
% inMaskPattern: input mask pattern; 1 for valid data and 0 for invalid data; 2D
% RegionRemoved: the size of region(RegionRemoved*RegionRemoved) to be removed; int
%
% outMaskpattern: output mask pattern; 1 for valid data and 0 for invalid data; 2D
% numRegionRemoved; number of removed region; can overlap; int 

% roughly acquiring biggest datas
data=sort(max(inPattern),'descend'); 
% roughly acquiring average biggest valid datas and setting the standard date
standard=mean(data(50:100))*3;
% acquirng the position of invalid data
[row,col]=find(inPattern>standard);

% the number of pixel; 1D
numPixel=size(inPattern,1);

for ii=1:size(row)
    rowRegion=max(1,row(ii)-RegionRemoved):min(row(ii)+RegionRemoved,numPixel);
    colRegion=max(1,col(ii)-RegionRemoved):min(col(ii)+RegionRemoved,numPixel);
    inMaskPattern(rowRegion,colRegion)=0;
end

outMaskpattern=inMaskPattern;
numRegionRemoved=size(row);
    
end