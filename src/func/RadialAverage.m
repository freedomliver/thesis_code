function [I,NumI,Std,PatternModified]=RadialAverage(pattern,correct)
% Radial average to each radial bin of 1-pixel width. 
% The center of the circle is the center of the pattern
% Any pixel that are more than (3) standard deviations away from the mean
% is removed
% 2020/11/5 qifengfeng
%
% pattern: m*n array of a single pattern
% correct: the factor of standard deviations away from the mean (usually set as 3)
%
% I: the radial average intense 
% NumI: the number of piexls of different radius
% Std: standard deviations for different radius
% PatternModified: any pixels that are more than (3) standard deviations
% away from the mean  are set as NaN

[y,x]=size(pattern); % the row index is y and the column index is x
% the center of the circle is the center of the pattern
centery=ceil(y/2);
centerx=ceil(x/2);

% Deciding the number of diferent radius
Rmax=max([sqrt((centery-1)^2+(centerx-1)^2),sqrt((centery-1)^2+(centerx-x)^2),...
    sqrt((centery-y)^2+(centerx-1)^2),sqrt((centery-y)^2+(centerx-x)^2)]);
Rmax=round(Rmax)+1;

I=zeros(1,Rmax);
NumI=I;
Std=I;
tmpI=I;
tmpNumI=I;
PatternModified=pattern;
Radius=zeros(y,x);

% Calculating the mean value for each radius
for i=1:y
    for j=1:x
        Radius(i,j)=round(sqrt((i-centery)^2+(j-centerx)^2)); 
        if Radius(i,j)==0 % +1: to avoid the apperance of index 0
            Radius(i,j)=1;
        end
        radius=Radius(i,j);
        deltaI= pattern(i,j);
        if ~isnan(deltaI) % remove the data of NaN
             tmpNumI(radius)=tmpNumI(radius)+1;
             tmpI(radius)=tmpI(radius)+deltaI;
        end
    end
end
tmpI=tmpI./tmpNumI;

% Calculating the standard deviation for each radius
for i=1:y
    for j=1:x
        radius=Radius(i,j);
        deltaI= pattern(i,j);
        if ~isnan(deltaI) % remove the data of NaN
              Std(radius)=Std(radius)+(deltaI-tmpI(radius))^2;
        end
    end
end
Std=sqrt(Std./tmpNumI);

% Any pixels that are more than (3) standard deviations away from the mean
for i=1:y
    for j=1:x
        radius=Radius(i,j);
        deltaI= pattern(i,j);
        % remove the data away from the mean more than (3) std;
        % 1 is reserved and 0 is set as NaN
        % the data of NaN is also included; the index of NaN is 0
        if (deltaI<=(tmpI(radius)+correct*Std(radius))) && (deltaI>=(tmpI(radius)-correct*Std(radius))) 
              NumI(radius)=NumI(radius)+1;    
              I(radius)=I(radius)+deltaI;
        else
              PatternModified(i,j)=NaN;
        end
    end
end
I=I./NumI;

Std=zeros(1,Rmax); % initialization
% Calculating the new standard deviation for each radius
for i=1:y
    for j=1:x
        radius=Radius(i,j);
        deltaI= PatternModified(i,j);
        if ~isnan(deltaI) % remove the data of NaN
              Std(radius)=Std(radius)+(deltaI-I(radius))^2;
        end
    end
end
Std=sqrt(Std./NumI);

end