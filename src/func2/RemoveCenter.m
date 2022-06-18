function [PatternOut] = RemoveCenter(PatternIn,radius)
%{
Remove the circle center of simulated diffraction pattern
2021/08/26 qifengfeng

PatternIn
radiu
PatternOut

%}

[row,col]=size(PatternIn);
center_row=round(row/2);
center_col=round(col/2);

PatternOut=PatternIn;
for ii=1:row
    for jj=1:col
        if norm([ii-center_row,jj-center_col]) <= radius
            PatternOut(ii,jj)=0;
        end
    end
end

end