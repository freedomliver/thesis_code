function psf_pattern = getPSFpattern(IM_Pattern, AngularDis, psf)
%GETPSFPATTERN��   IM_Patternͨ��PSF������ӽǷֲ�AngularDis�µ�����ͼ��
%   �˴���ʾ��ϸ˵��
%  IM_Pattern����ʼ����ͼ
%  AngularDis�����ӽǷֲ�
%  psf��PSF
%  pat_size�� ����ͼ��С
%  psf_pattern������õ���������ͼ

pat_size = 400;  x = 1 : pat_size;  y = 1 : pat_size;   % ����ͼ��С
nub_theta = 41; nub_beta = 30;                         % na �ǦȽǡ�0���С���nb�Ǧսǡ�0,2�С�
theta = pi/nub_theta/2 : (pi/nub_theta) : pi-pi/nub_theta/2; 
[xx, yy] = meshgrid(-0.025*2*199-0.025 : 0.025*2:0.025*2*200-0.025, -0.025*2*199-0.025 : 0.025*2:0.025*2*200-0.025);
[~, RR] = cart2pol(xx, yy);       %��ֱ�������ɼ�����  
Dist_total=0;     %�����һ������
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