% % This script is prepared for thesis 2-3
% clc;  clear;  close all;
% addpath(genpath('D:\Program Files\MATLAB\cal\src'));    % Add the path of program package
% DPWAfile = 'D:\Program Files\MATLAB\cal\src\DPWA\';   % Diffraction amplitude path (DPWA file)
% sx = 0 : 0.025 : 10;  sy = -10 : 0.025 : 10;
% smap = zeros(length(sy), length(sx));      % scattering amplitude
% for ii = 1 : length(sy), for jj = 1 : length(sx), smap(ii, jj) = norm([sy(ii), sx(jj)]); end; end
% [fH, ~] = DPWAcpu(smap, 1, 3000, DPWAfile);
% [fHe,~]= DPWAcpu(smap, 2, 3000, DPWAfile);
% [fC, ~] = DPWAcpu(smap, 6, 3000, DPWAfile);
% [fN, ~] = DPWAcpu(smap, 7, 3000, DPWAfile);
% [fO, ~] = DPWAcpu(smap, 8, 3000, DPWAfile);
% [fF, ~ ] = DPWAcpu(smap, 9, 3000, DPWAfile);
% [fS, ~] = DPWAcpu(smap, 16,3000, DPWAfile);
% [fAr, ~]=DPWAcpu(smap, 18,3000, DPWAfile);
% [fCr, ~]= DPWAcpu(smap, 24,3000, DPWAfile);
% [fI, ~ ] = DPWAcpu(smap, 53,3000, DPWAfile);
% 
% %  N2  
% rN1 = [0; 0.5488; 0];     % Aligned along y axis
% rN2 = [0; -0.5488; 0];
% molecular = [rN1, rN2];

%% Euler angle
cossquare = load('D:\Program Files\MATLAB\cal\thesis_code\Cossquare-N2-Engery_density_1.569-Temp_30K.dat');
Evolution = load('D:\Program Files\MATLAB\cal\thesis_code\Evolution-N2-Engery_density_1.569-Temp_30K.dat');
[m1, p1] = max(cossquare);  [m2, p2] = min(cossquare);
distri_max = Evolution(p1(2), :); distri_min = Evolution(p2(2), :);
beta = linspace(0, pi, length(distri_max));

alpha = 0 : (2*pi/60) : 2*pi;     
gamma = 0;
% AngularDis_min = distri_min .* sin(beta);       
% AngularDis_max = distri_max .* sin(beta);       
% AngularDis_rand = sin(beta);

% figure, plot(beta, AngularDis, '-x'); grid on
% xlabel('$\theta$', 'Interpreter', 'latex'),  ylabel('f($\theta$) (arb.unit)', 'Interpreter', 'latex');
% set(gca, 'linewidth', 1),  set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
% set(gca, 'Xtick', [0,pi/2,pi]),  set(gca, 'XTicklabel', {'0','\pi/2','\pi'});

%%
IT = smap .* 0;  IA = smap .* 0; 
for beta_index = 1 : length(beta) 
    theta = beta(beta_index);
    layer1 = repmat(sx, length(sy), 1);       % sx
    layer2 = repmat(sy', 1, length(sx));      % sy
    layer3 = zeros(length(sy), length(sx)); % sz=0
    s = cat(3, layer1, layer2, layer3);          % dimension3: [sx,sy,sz]
    
    It = smap .* 0;  Ia = smap .* 0;
    for alpha_index = 1 : length(alpha)
        phi = alpha(alpha_index);       
        for gamma_index = 1 : length(gamma) 
            roll = gamma(gamma_index);
            r = RotationMatrixYZY(phi, theta, roll) * molecular;          
            Ia = Ia + abs(fN).^2.*2;      
            diffraciton = fN .* exp(1i.*sum(s.*reshape(r(:,1), [1,1,3]), 3)) + ...
                                  fN .* exp(1i.*sum(s.*reshape(r(:,2), [1,1,3]), 3));   
            It = It + abs(diffraciton).^2;             
        end
    end
%     AngularDis_max = 1;
    IT = IT + It ./ (length(alpha) .* length(gamma)) .* AngularDis_min(beta_index);
    IA = IA+ Ia ./ (length(alpha) .* length(gamma)) .* AngularDis_min(beta_index);
end
IA = IA ./ sum(AngularDis_min);
IT = IT ./ sum(AngularDis_min);

IT_Pattern = [fliplr(IT(:, 2:end)), IT] ./ max(IT, [], 'all');
IA_Pattern = [fliplr(IA(:, 2:end)), IA] ./ max(IA, [], 'all');
IM = IT - IA;  IM_Pattern = [fliplr(IM(:, 2:end)), IM] ./ max(IT, [], 'all');
M = IM ./ IA;  M_left = fliplr(M(:, 2:end));  M_Pattern = [M_left, M];

figure
pcolor(sy,sy,M_Pattern), colorbar; shading interp,colormap jet
xlabel('sx ($\AA^{-1}$)','Interpreter','latex'), ylabel('sy ($\AA^{-1}$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',20)
set(gca,'XTick',-12:4:12),set(gca,'XTicklabel',-12:4:12)
set(gca,'YTick',-12:4:12),set(gca,'YTicklabel',-12:4:12)

figure;
imagesc(sy, sy, IT_Pattern);  colorbar;
xlabel('sx ($\AA^{-1}$)','Interpreter','latex'), ylabel('sy ($\AA^{-1}$)','Interpreter','latex')
set(gca, 'FontSize', 25);  set(gcf, 'unit', 'centimeters', 'position', [5 2 18 18]);  set(gca, 'Position', [.15 .2 .65 .65]);
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 

IT_Pattern_min = IT_Pattern; 
IA_Pattern_min = IA_Pattern;
M_Pattern_min = M_Pattern;
% 
% IT_Pattern_rand = IT_Pattern;
% IA_Pattern_rand = IA_Pattern;
% M_Pattern_rand = M_Pattern;

% IT_Pattern_max = IT_Pattern;
% IA_Pattern_max = IA_Pattern;
% M_Pattern_max = M_Pattern;

%% Fourier and Hankel transf of N2 with Perfect alignment and Random Orientation
attenu_pat_min = M_Pattern_min;  attenu_pat_rand = M_Pattern_rand; attenu_pat_max = M_Pattern_max; 
pattern_size = size(attenu_pat_min); 
r_size = (pattern_size(1)-1) / 2;  delta_s = 0.05;
for mm = 1 : pattern_size(1)
    for nn = 1 : pattern_size(1)
        s = sqrt((mm-r_size-1)^2 + (nn-r_size-1)^2) * delta_s;
        attenu_pat_min(mm, nn) = attenu_pat_min(mm, nn) * exp(-1*s^2 / 6.5^2);
        attenu_pat_rand(mm, nn) = attenu_pat_rand(mm, nn) * exp(-1*s^2 / 6.5^2);       
        attenu_pat_max(mm, nn) = attenu_pat_max(mm, nn) * exp(-1*s^2 / 6.5^2);
    end
end

signal_len = 2^13;         % length of signal, pad with 0
delta_k = delta_s / (2*pi);    % sampling period 
sampling_fre = 1 / delta_k;  % sampling frequency
delta_r = sampling_fre / signal_len;

rand_pat = zeros(signal_len, signal_len);  min_pat = rand_pat;  max_pat = rand_pat;
min_pat((signal_len/2+1-r_size) : (signal_len/2+1+r_size),...
              (signal_len/2+1-r_size) : (signal_len/2+1+r_size)) = attenu_pat_min;
rand_pat((signal_len/2+1-r_size) : (signal_len/2+1+r_size),...
              (signal_len/2+1-r_size) : (signal_len/2+1+r_size)) = attenu_pat_rand;
max_pat((signal_len/2+1-r_size) : (signal_len/2+1+r_size),...
              (signal_len/2+1-r_size) : (signal_len/2+1+r_size)) = attenu_pat_max;

pattern_len = 300;
pdf_min = FourierHankel(min_pat, signal_len, pattern_len);
pdf_rand = FourierHankel(rand_pat, signal_len, pattern_len);
pdf_max = FourierHankel(max_pat, signal_len, pattern_len);
pdf_min = pdf_min ./ max(pdf_min,[],'all');
pdf_rand = pdf_rand ./ max(pdf_rand,[],'all');
pdf_max = pdf_max ./ max(pdf_max,[],'all');
rx = (-1*pattern_len : pattern_len) .* delta_r;

save('article_chapter2.3.1.mat');

%% figure 2.3 IM, sM and fr picture for different N2 alignment
load('2.3_N2_aligned_diffraction.mat');
fig = figure; set(gcf, 'position', [700, 90, 1000, 900]);    % subplot[3,3]. IM, sM and f(r) for each side line
ax(11) = axes('Position', [0.08, 0.69, 0.4, 0.264]);
ax(12) = axes('Position', [0.08, 0.37, 0.4, 0.262]);
ax(13) = axes('Position', [0.08, 0.055, 0.4, 0.26]);
% ax(21) = axes('Position', [0.39, 0.69, 0.27, 0.264]);
% ax(22) = axes('Position', [0.39, 0.37, 0.27, 0.262]);
% ax(23) = axes('Position', [0.39, 0.055, 0.27, 0.26]);
ax(31) = axes('Position', [0.55, 0.69, 0.4, 0.264]);
ax(32) = axes('Position', [0.55, 0.37, 0.4, 0.262]);
ax(33) = axes('Position', [0.55, 0.055, 0.4, 0.26]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

imagesc(ax(11), sy, sy, IT_Pattern_min);  colorbar(ax(11));text(ax(11), -9, -9, '(A)', 'color', 'w');
xlabel(ax(11), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'), ylabel(ax(11), 'sy ($\AA^{-1}$)', 'Interpreter', 'latex');
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
title(ax(11), 'cos^2θ = 0.2044对应的IT、sM和PDF');
set(ax(11), 'FontSize', 12);

imagesc(ax(12), sy, sy, M_Pattern_min);  colorbar(ax(12));text(ax(12), -9, -9, '(B)', 'color', 'w');
xlabel(ax(12), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'), ylabel(ax(12), 'sy ($\AA^{-1}$)', 'Interpreter', 'latex');
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(12), 'FontSize', 12);

imagesc(ax(13), rx, rx, pdf_min), colorbar(ax(13)); text(ax(13), -4, -4, '(C)', 'color', 'w');
xlabel(ax(13), 'rx ($\AA$)', 'Interpreter', 'latex'), ylabel(ax(13), 'ry ($\AA$)', 'Interpreter', 'latex');
set(gca, 'xtick', -3 : 3), set(gca, 'ytick', -3 : 3);
set(ax(13), 'FontSize', 12);

% imagesc(ax(21), sy, sy, IT_Pattern_rand);  colorbar(ax(21));
% xlabel(ax(21), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'), %ylabel(ax(21), 'sy ($\AA^{-1}$)', 'Interpreter', 'latex');
% set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
% set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
% title(ax(21), 'cos^2θ = 0.33对应的IT、sM和PDF');
% set(ax(21), 'FontSize', 12);
% 
% imagesc(ax(22), sy, sy, M_Pattern_rand);  colorbar(ax(22));
% xlabel(ax(22), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'),% ylabel(ax(22), 'sy ($\AA^{-1}$)', 'Interpreter', 'latex');
% set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
% set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
% set(ax(22), 'FontSize', 12);
% 
% imagesc(ax(23), rx, rx, pdf_rand), colorbar(ax(23)); 
% xlabel(ax(23), 'rx ($\AA$)', 'Interpreter', 'latex'), ylabel(ax(23), 'ry ($\AA$)', 'Interpreter', 'latex');
% set(gca, 'xtick', -3 : 3), set(gca, 'ytick', -3 : 3);
% set(ax(23), 'FontSize', 12);

imagesc(ax(31), sy, sy, IT_Pattern_max);  colorbar(ax(31));  text(ax(31), -9, -9, '(a)', 'color', 'w');
xlabel(ax(31), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'), %ylabel(ax(31), 'sy ($\AA^{-1}$)', 'Interpreter', 'latex');
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
title(ax(31), 'cos^2θ = 0.4255对应的IT、sM和PDF');
set(ax(31), 'FontSize', 12);

imagesc(ax(32), sy, sy, M_Pattern_max);  colorbar(ax(32)); text(ax(32), -9, -9, '(b)', 'color', 'w');
xlabel(ax(32), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'),% ylabel(ax(32), 'sy ($\AA^{-1}$)', 'Interpreter', 'latex');
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(32), 'FontSize', 12);

imagesc(ax(33), rx, rx, pdf_max), colorbar(ax(33));   text(ax(33), -4, -4, '(c)', 'color', 'w');
xlabel(ax(33), 'rx ($\AA$)', 'Interpreter', 'latex'), ylabel(ax(33), 'ry ($\AA$)', 'Interpreter', 'latex');
set(gca, 'xtick', -3 : 3), set(gca, 'ytick', -3 : 3);
set(ax(33), 'FontSize', 12);

%%
cossquare = load('D:\Program Files\MATLAB\cal\thesis_code\Cossquare-N2-Engery_density_1.569-Temp_30K.dat');
Evolution = load('D:\Program Files\MATLAB\cal\thesis_code\Evolution-N2-Engery_density_1.569-Temp_30K.dat');
figure;
plot(cossquare, '-x');
% [m1, p1] = max(cossquare);  [m2, p2] = min(cossquare);
% distri_max = Evolution(p1(2), :); distri_min = Evolution(p2(2), :);
% beta = linspace(0, pi, length(distri_max));
%计算完美准直的信号，需要重新计算，按照荧光屏的转移动量大小计算
%{
beta = 0;
alpha = 0 : (2*pi/60) : 2*pi;     
gamma = 0;
IT = smap .* 0;  IA = smap .* 0; 
for beta_index = 1 : length(beta) 
    theta = beta(beta_index);
    layer1 = repmat(sx, length(sy), 1);       % sx
    layer2 = repmat(sy', 1, length(sx));      % sy
    layer3 = zeros(length(sy), length(sx)); % sz=0
    s = cat(3, layer1, layer2, layer3);          % dimension3: [sx,sy,sz]
    
    It = smap .* 0;  Ia = smap .* 0;
    for alpha_index = 1 : length(alpha)
        phi = alpha(alpha_index);       
        for gamma_index = 1 : length(gamma) 
            roll = gamma(gamma_index);
            r = RotationMatrixYZY(phi, theta, roll) * molecular;          
            Ia = Ia + abs(fN).^2.*2;      
            diffraciton = fN .* exp(1i.*sum(s.*reshape(r(:,1), [1,1,3]), 3)) + ...
                                  fN .* exp(1i.*sum(s.*reshape(r(:,2), [1,1,3]), 3));   
            It = It + abs(diffraciton).^2;             
        end
    end
    IT = IT + It ./ (length(alpha) .* length(gamma));
    IA = IA+ Ia ./ (length(alpha) .* length(gamma));
end
IT_Pattern = [fliplr(IT(:, 2:end)), IT] ./ max(IT, [], 'all');
IA_Pattern = [fliplr(IA(:, 2:end)), IA] ./ max(IA, [], 'all');
IM = IT - IA;  IM_Pattern = [fliplr(IM(:, 2:end)), IM] ./ max(IT, [], 'all');
M = IM ./ IA;  M_left = fliplr(M(:, 2:end));  M_Pattern = [M_left, M];

figure
pcolor(sy,sy,M_Pattern), colorbar; shading interp,colormap jet
xlabel('sx ($\AA^{-1}$)','Interpreter','latex'), ylabel('sy ($\AA^{-1}$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',20)
set(gca,'XTick',-12:4:12),set(gca,'XTicklabel',-12:4:12)
set(gca,'YTick',-12:4:12),set(gca,'YTicklabel',-12:4:12)
%}   
%% 计算PSF

pat_size = 400;  x = 1 : pat_size;  y = 1 : pat_size;   % 衍射图大小
nub_theta = 41; nub_beta = 30;                         % na 是θ角【0，π】，nb是φ角【0,2π】
theta = pi/nub_theta/2 : (pi/nub_theta) : pi-pi/nub_theta/2; 
beta = 2*pi/nub_beta/2 : (2*pi/nub_beta) : 2*pi-2*pi/nub_beta/2;     
psf = cell(pat_size, pat_size);                            % psf=zeros(n,n,na,nb,2);  

[xx, yy] = meshgrid(-0.025*2*199-0.025 : 0.025*2:0.025*2*200-0.025, -0.025*2*199-0.025 : 0.025*2:0.025*2*200-0.025);
[Phib, RR] = cart2pol(xx, yy);       %把直角坐标变成极坐标  
% Phib = rot90(Phib, 3);           % 为了
tic;
for index_theta = 1 : nub_theta
    for index_beta = 1 : nub_beta
        for i = 1 : pat_size
            for j = 1 : pat_size
                if RR(j,i) < pat_size/2 .* 0.025*2         % 最大的圆的范围内
                    phib = Phib(j,i);                    %  两个角度的值, 探测网格上的角度
                    phi = acos(sin(theta(index_theta)) .* cos(beta(index_beta)) .* sin(phib) + cos(theta(index_theta)) .* cos(phib));
                     %       acos(sinθ .* cosφ .* sin(meshΘb)) + cosθ .* cos(meshΘb)
                    psfx = round(RR(j,i) .* sin(phi) ./ (0.025*2) + length(x)/2 + 0.5);
                    psfy = round(RR(j,i) .* cos(phi) ./ (0.025*2) + length(y)/2 + 0.5);
                    psf{j,i}(index_theta, index_beta, 2) = psfx; 
                    psf{j,i}(index_theta, index_beta, 1) = psfy;
                end
            end
        end
    end
    toc
end

save('PSF_patsize400_theta41.mat');

%% 测试 对于完美准直衍射信号，沿着几乎完美准直角分布PSF  通过PSF计算得到的衍射信号
load('PSF_patsize400_theta41.mat');
time_some = 1:10:3550;
GET_N2_pattern_diff = zeros(pat_size, pat_size,length(time_some));  tic
% IM_PATTERN_N2_rand = getPSFpattern(IM_Pattern, ones(1,41), psf);
for i = 1 : length(time_some)
    angulaDis_t = time_some(i); 
    GET_N2_pattern_diff(:,:, i)=getPSFpattern(IM_Pattern, Evolution(angulaDis_t,:), psf) - IM_PATTERN_N2_rand; 
    toc 
end

%%  计算各项异性值的部分
cc = 201;   bias = 101;  a1 = 59;  b1 = 31;    
x1 = cc-bias-a1;  x2 = cc-bias+a1; x3 = cc+bias-a1;  x4 = cc+bias+a1;   y1 = cc-b1; y2 = cc+b1;

backsize = size(IM_PATTERN_N2_rand);
diff0 = zeros(1, length(time_some));
for i=1:length(time_some)
    figure(97);  imagesc(GET_N2_pattern_diff(:, :, i)); title(['i=',num2str(i)]); pause(0.001);   colorbar;  hold on 
    plot([ x1,x1,x2,x2,x3,x3,x4,x4,y1,y2,y1,y2,y1,y2,y1,y2],...
            [y1,y2,y1,y2,y1,y2,y1,y2,x1,x1,x2,x2,x3,x3,x4,x4], 'xr', 'MarkerSize', 15);
    
    ave_a = mean(GET_N2_pattern_diff(x1:x2, y1:y2, i), 'all');
    ave_b = mean(GET_N2_pattern_diff(x3:x4, y1:y2, i),'all');
    ave_c = mean(GET_N2_pattern_diff(y1:y2, x1:x2, i),'all');
    ave_d = mean(GET_N2_pattern_diff(y1:y2, x3:x4, i), 'all');
    diff0(i) = ave_a + ave_b - ave_c - ave_d;
end

figure;
plot(time_some, diff0);

%%  氮气各向异性模拟图
load('2.3_N2_anistropy.mat');
c = 2.9979e8;  hbar = 1.054589e-34;  B0 = 1.998;  B = B0 * hbar * 2*pi * c * 1e2; 
time = time_some * hbar / B * 1e9;
align_t = cossquare(:,1);  align_cossquare = cossquare(:,2);
figure;
yyaxis left; plot(time, diff0); ylabel('anisotropy'); hold on; yyaxis right; 
plot(align_t(1:10:3550),align_cossquare(1:10:3550),'-x'); ylabel('<cos^2\theta>');
xlabel('times(ps)');title('氮气分子准直度和各向异性模拟图'); set(gca, 'FontSize', 13);
axes('position', [0.20, 0.65, 0.18,0.25]);
imagesc(GET_N2_pattern_diff(:, :, 300));  hold on
plot([ x1,x1,x2,x2,x3,x3,x4,x4,y1,y2,y1,y2,y1,y2,y1,y2],...
            [y1,y2,y1,y2,y1,y2,y1,y2,x1,x1,x2,x2,x3,x3,x4,x4], 'xr', 'MarkerSize', 4);
set(gca, 'xtick', [], 'xticklabel', []);set(gca, 'ytick', [], 'yticklabel', []);

