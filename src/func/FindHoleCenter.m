function [HoleCenterRow,HoleCenterCol,HoleRadius] = FindHoleCenter(pattern,HoleStrengthMax,HoleDeviation)
% find the center of the hole in fluorescent screen; all patterns have the same hole;
% set the mask pattern in this region as 0
% Here can be replaced by circlar hough transform
% 2020/10/15 qifengfeng
%
% pattern: the pattern prepared for finding the hole
% HoleStrengthMax: the upper limit of the pattern strength in the hole; for converting the pattern to a binary image
% HoleDeviation: the rough pixels of the hole edge off the center of the pattern
%
% HoleCenterRow: the row center of the hole(y);int
% HoleCenterCol: the col center of the hole(x); int
% HoleRadius: the radius of the hole

[numRowPix,numColPix]=size(pattern); % acquire the size of the pattern; y for row and x for column

% converting the pattern to a binary image
pattern=pattern<HoleStrengthMax;
% check the pattern of the hole
RowData=(round(numRowPix/2)-HoleDeviation):(round(numRowPix/2)+HoleDeviation);
ColData=(round(numColPix/2)-HoleDeviation):(round(numColPix/2)+HoleDeviation);

subplot(1,2,1)
imshow(pattern(RowData,ColData))

% acquring the center and radius
hole = regionprops(pattern(RowData,ColData),'Centroid','MajorAxisLength','MinorAxisLength');
% all radius 
radius=cat(1,hole.MajorAxisLength); 
% the accurate hole
[~,pos]=max(radius);

HoleCenterCol=round(round(numColPix/2)-HoleDeviation+hole(pos).Centroid(1)-1); % X; col
HoleCenterRow=round(round(numRowPix/2)-HoleDeviation+hole(pos).Centroid(2)-1); % Y; row
HoleRadius=round((hole(pos).MajorAxisLength+hole(pos).MinorAxisLength)/2/2);

subplot(1,2,2)
imshow(pattern)
hold on 
viscircles([HoleCenterCol,HoleCenterRow],HoleRadius,'EdgeColor','r');
set(gcf,'unit','centimeters','position',[5,10,30,15])

disp('===the hole in fluorescent screen has been found===')

end