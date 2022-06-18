function [polor_pattern,xyz_pattern] = Cartesian2polor_recover(p_r,Legendre_order,max_r,theta,radius)
% converting pattern form Cartesian coordinates to polar coordinates
% Legendre fitting to order 10
% recover the pattern using Legendre fitting result
% used with funciton 'Cartesian2polor_fit'
%
% p_r: fitting result
% Legendre_order: reserved Legendre_order
% max_r: limited radius
% theta: theta range of the pattern
% radius: radius range of the pattern
%
% polor_pattern: pattern in polor coordinate
% xyz_pattern: pattern in xyz coordinate

Legendre_num=length(Legendre_order);

% xyz_pattern
xyz_pattern=zeros(2*max_r+1,2*max_r+1);
for ii=1:2*max_r+1
    for jj=1:2*max_r+1
        rr=round(((ii-max_r-1).^2+(jj-max_r-1).^2).^0.5);
        cos_theta=(max_r+1-ii)/rr;
        if rr>0 && rr<=max_r
            for index=1:Legendre_num
                xyz_pattern(ii,jj)=xyz_pattern(ii,jj)+...
                                   p_r(Legendre_order(index)+1,rr)*LegendreZero(Legendre_order(index),cos_theta);
            end
        end
    end
end

% polor_pattern
num_theta=length(theta);
num_radius=length(radius);
polor_pattern=zeros(num_theta,num_radius);

for ii=1:num_theta
    for jj=1:num_radius
        for index=1:Legendre_num
            polor_pattern(ii,jj)=polor_pattern(ii,jj)+...
                 p_r(Legendre_order(index)+1,radius(jj))*LegendreZero(Legendre_order(index),cos(theta(ii)));
        end
    end
end

end
