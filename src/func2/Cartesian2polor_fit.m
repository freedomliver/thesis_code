function [p_r] = Cartesian2polor_fit(pattern,max_r)
% converting pattern form Cartesian coordinates to polar coordinates
% Legendre fitting to order 10
% used with funciton 'Cartesian2polor_recover'
%
% pattern: pattern to be processed
% max_r: limited radius
%
% p_r: fitting result

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

end
