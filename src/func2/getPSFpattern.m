function psf_pattern = getPSFpattern(IM_Pattern, AngularDis, psf)
%GETPSFPATTERN：   IM_Pattern通过PSF计算分子角分布AngularDis下的衍射图像
%   此处显示详细说明
%  IM_Pattern：初始衍射图
%  AngularDis：分子角分布
%  psf：PSF
%  pat_size： 衍射图大小
%  psf_pattern：计算得到的新衍射图

pat_size = 400;  x = 1 : pat_size;  y = 1 : pat_size;   % 衍射图大小
nub_theta = 41; nub_beta = 30;                         % na 是θ角【0，π】，nb是φ角【0,2π】
theta = pi/nub_theta/2 : (pi/nub_theta) : pi-pi/nub_theta/2; 
[xx, yy] = meshgrid(-0.025*2*199-0.025 : 0.025*2:0.025*2*200-0.025, -0.025*2*199-0.025 : 0.025*2:0.025*2*200-0.025);
[~, RR] = cart2pol(xx, yy);       %把直角坐标变成极坐标  
Dist_total=0;     %计算归一化因子
normfac = zeros(1, nub_theta);
for index_theta = 1 : nub_theta       
    normfac(index_theta) = sin(theta(index_theta)) .* AngularDis(index_theta);
    Dist_total = Dist_total + sin(theta(index_theta)) .* AngularDis(index_theta);
end

Mperf0 = IM_Pattern;
Orientation = psf{1,1};
psf_pattern = zeros(pat_size);
for i = 1 : pat_size    %n
    for j = 1 : pat_size
        Orientation = psf{j,i};
        if RR(j,i) < pat_size/2 .* 0.025*2 & RR(j,i) > pat_size/20/2 .* 0.025*2
            for index_theta = 1 : nub_theta
                for index_beta = 1 : nub_beta
                    psf_pattern(j, i) = psf_pattern(j, i) + Mperf0(Orientation(index_theta, index_beta, 1), ...
                                                                                             Orientation(index_theta, index_beta, 2)) .* normfac(index_theta);
                end
            end
        end
    end
end
psf_pattern = psf_pattern ./ nub_beta ./ Dist_total;

end