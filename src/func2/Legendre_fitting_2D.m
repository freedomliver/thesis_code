function [pattern_out] = Legendre_fitting_2D(pattern,max_r,Legendre_order)
% 2D Legendre fitting
%
% pattern: pattern to be processed
% max_r: limited radius of the pattern
% Legendre_order: reserved Legendre order when recovering the pattern
%
% pattern_out: pattern output

% Legendre fitting
[row,col]=size(pattern);
radius=zeros(row,col);
cos_theta=zeros(row,col);
center=[round(row/2),round(col/2)];

pattern_polor=cell(max_r,1);
for ii=1:row
    for jj=1:col
        radius(ii,jj)=round(((ii-center(1)).^2+(jj-center(2)).^2).^0.5);
        cos_theta(ii,jj)=(center(1)-ii)/radius(ii,jj);
        if radius(ii,jj)>0 && radius(ii,jj)<=max_r
        pattern_polor{radius(ii,jj),1}(end+1,:)=[cos_theta(ii,jj),pattern(ii,jj)];
        end
    end
end

p_r=zeros(10,max_r);
for ii=10:max_r
    [p,~]=Legendre_fit(pattern_polor{ii,1}(:,1),pattern_polor{ii,1}(:,2));
    p_r(1,ii)=p.a; p_r(2,ii)=p.b;p_r(3,ii)=p.c;p_r(4,ii)=p.d;p_r(5,ii)=p.e;
    p_r(6,ii)=p.f;p_r(7,ii)=p.g;p_r(8,ii)=p.h;p_r(9,ii)=p.m;p_r(10,ii)=p.n;
end

% recovering the pattern
Legendre_num=length(Legendre_order);
pattern_out=zeros(2*max_r+1,2*max_r+1);
for ii=1:2*max_r+1
    for jj=1:2*max_r+1
        rr=round(((ii-max_r-1).^2+(jj-max_r-1).^2).^0.5);
        cos_theta=(max_r+1-ii)/rr;
        if rr>0 && rr<=max_r
            for index=1:Legendre_num
                pattern_out(ii,jj)=pattern_out(ii,jj)+...
                    p_r(Legendre_order(index)+1,rr)*LegendreZero(Legendre_order(index),cos_theta);
            end
        end
    end
end

end
