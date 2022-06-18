function [MeanRegion,FitRegion]=PreFindPatternCenter(FindCenterRadius,RadiusCoefficient,SizePix,HoleCenter)
% creating the region for calculating the mean data as standard for circle fitting
% and creating the region for fitting the circle
% pre information for the function FindPatternCenter
% avoiding repetitive computation
% 2020/10/20 qifengfeng
%
% FindCenterRadius: the radius to be fitted; 1-D array
% RadiusCoefficient; the coefficient of the radius to confirm the range of the fitted region; such as [0.9,1.1]
% SizePix: the size of the pattern; [row number, col number]
% HoleCenter: the center of the hole in the pattern; rough center of the pattern; [row,col]
%
% MeanRegion: the region for calculating the mean data as standard
% FitRegion: the region for fitting the circle

% number of the radius to be fitted
numRadius=length(FindCenterRadius);
% the region for calculating the mean data as standard
MeanRegion=zeros(SizePix(1),SizePix(2),numRadius); 
% the region for fitting the circle
FitRegion=zeros(SizePix(1),SizePix(2),numRadius);

for ii=1:numRadius
    Outer=CircleRegion(SizePix,FindCenterRadius(ii)+3,HoleCenter);
    Inner=CircleRegion(SizePix,FindCenterRadius(ii)-3,HoleCenter);
    MeanRegion(:,:,ii)=Outer.*~Inner;
    
    Outer=CircleRegion(SizePix,FindCenterRadius(ii)*RadiusCoefficient(2),HoleCenter);
    Inner=CircleRegion(SizePix,FindCenterRadius(ii)*RadiusCoefficient(1),HoleCenter);
    FitRegion(:,:,ii)=Outer.*~Inner;
end

end

function Circle = CircleRegion(SizePix,Radius,HoleCenter)
% setting datas in a circle region of matrix as 1, others are 0
%
% SizePix: the size of the pattern; [row number, col number]
% Radius: the radius of the circle
% HoleCenter: the center of the hole in the pattern; rough center of the pattern; [row,col]
%
% Circle: the matrix where datas  in a circle region are 1 and others are 0

Circle=zeros(SizePix); % initial the matrix
for Row=max(1,HoleCenter(1)-Radius):min(SizePix(1),HoleCenter(1)+Radius)
    for Col=max(1,HoleCenter(2)-Radius):min(SizePix(2),HoleCenter(2)+Radius)
        if ((Row-HoleCenter(1))^2+(Col-HoleCenter(2))^2)<=(Radius^2)
            Circle(Row,Col)=1;
        end
    end
end

end