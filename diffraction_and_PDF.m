% This script is prepare for my undergraduate thesis second chapter
addpath(genpath('D:\Program Files\MATLAB\cal\src'));    % Add the path of program package
DPWAfile = 'D:\Program Files\MATLAB\cal\src\DPWA\';   % Diffraction amplitude path (DPWA file)
sx = -10 : 0.05 : 10;  sy = -10 : 0.05 : 10;                                 % Changable according to used data-panel
smap = zeros(length(sy), length(sx));      % scattering amplitude
for ii = 1 : length(sy), for jj = 1 : length(sx), smap(ii, jj) = norm([sy(ii), sx(jj)]); end; end
% [fH, ~] = DPWAcpu(smap, 1, 3000, DPWAfile);
% [fHe,~]= DPWAcpu(smap, 2, 3000, DPWAfile);
% [fC, ~] = DPWAcpu(smap, 6, 3000, DPWAfile);
[fN, ~] = DPWAcpu(smap, 7, 3000, DPWAfile);
% [fO, ~] = DPWAcpu(smap, 8, 3000, DPWAfile);
% [fF, ~ ] = DPWAcpu(smap, 9, 3000, DPWAfile);
% [fS, ~] = DPWAcpu(smap, 16,3000, DPWAfile);
% [fAr, ~]=DPWAcpu(smap, 18,3000, DPWAfile);
% [fCr, ~]= DPWAcpu(smap, 24,3000, DPWAfile);
% [fI, ~ ] = DPWAcpu(smap, 53,3000, DPWAfile);

%% Chapter2.  Part 2.1, calculate the perfert molecule differaction for N2
% N2 molecular information
rN1 = [0; 0.5488; 0];     % Aligned along y axis
rN2 = [0; -0.5488; 0];
molecular = [rN1, rN2];

%% Perfect Alignment of N2  only one molecule
alpha = 0 : (2*pi/60) : 2*pi;     beta = 0;  gamma = 0;  % Euler angle
AngularDis = 1;   % angular distribution
 
% Calculate the diffraction
    layer1 = repmat(sx, length(sy), 1);       % sx
    layer2 = repmat(sy', 1, length(sx));      % sy
    layer3 = zeros(length(sy), length(sx)); % sz=0
    s = cat(3, layer1, layer2, layer3);          % dimension3: [sx,sy,sz]
IT = smap .* 0;  IA = smap .* 0; 
    for alpha_index = 1 : length(alpha) % alpha
        phi = alpha(alpha_index);       
            r = RotationMatrixYZY(phi, beta, gamma) * molecular; % r             
            IA = IA + abs(fN).^2.*2;      
            diffraciton = fN .* exp(1i.*sum(s.*reshape(r(:,1), [1,1,3]), 3)) + ...
                                  fN .* exp(1i.*sum(s.*reshape(r(:,2), [1,1,3]), 3));   
            IT = IT + abs(diffraciton).^2;             
    end
IT = IT ./ sum(AngularDis) ./ length(alpha);  IT_Pattern_perf = IT ./ max(IT, [], 'all');
IA = IA ./ sum(AngularDis) ./ length(alpha); IA_Pattern_perf = IA ./ max(IA, [], 'all');
IM = IT - IA;   IM_Pattern_perf = IM ./ max(IT, [], 'all');
M_Pattern_perf = IM ./ IA;

%%   random orientation of N2   length(beta) * length(alpha) * length(gamma) molecules
alpha = 0 : (2*pi/60) : 2*pi;     beta = 0 : (pi/20) : pi;  gamma = 0;  % Euler angle
AngularDis = sin(beta);   % angular distribution
% Calculate the differaction
    layer1 = repmat(sx, length(sy), 1);       % sx
    layer2 = repmat(sy', 1, length(sx));      % sy
    layer3 = zeros(length(sy), length(sx)); % sz=0
    s = cat(3, layer1, layer2, layer3);          % dimension3: [sx,sy,sz]
IT = smap .* 0;  IA = smap .* 0; 
for beta_index = 1 : length(beta)
    ia = 0; it = 0;
    for alpha_index = 1 : length(alpha) % alpha
        phi = alpha(alpha_index);       
        % if gamma not = 0, put for gamma_index in here
            r = RotationMatrixYZY(phi, beta(beta_index), gamma) * molecular; % r             
            ia = ia + abs(fN).^2.*2;      
            diffraciton = fN .* exp(1i.*sum(s.*reshape(r(:,1), [1,1,3]), 3)) + ...
                                  fN .* exp(1i.*sum(s.*reshape(r(:,2), [1,1,3]), 3));   
            it = it + abs(diffraciton).^2;             
    end
    IA = IA + ia .* AngularDis(beta_index);
    IT = IT + it .* AngularDis(beta_index);
end
IT = IT ./ sum(AngularDis) ./ length(alpha);  IT_Pattern_rand = IT ./ max(IT, [], 'all');
IA = IA ./ sum(AngularDis) ./ length(alpha); IA_Pattern_rand = IA ./ max(IA, [], 'all');
IM = IT - IA;   IM_Pattern_rand = IM ./ max(IT, [], 'all');
M_Pattern_rand = IM ./ IA;

%% Fourier and Hankel transf of N2 with Perfect alignment and Random Orientation
% just the same with inves-abel transf
attenu_pat_perf = M_Pattern_perf;  attenu_pat_rand = M_Pattern_rand; 
pattern_size = size(attenu_pat_perf); 
r_size = (pattern_size(1)-1) / 2;  delta_s = 0.05;
for mm = 1 : pattern_size(1)
    for nn = 1 : pattern_size(1)
        s = sqrt((mm-r_size-1)^2 + (nn-r_size-1)^2) * delta_s;
        attenu_pat_perf(mm, nn) = attenu_pat_perf(mm, nn) * exp(-1*s^2 / 6.5^2);
        attenu_pat_rand(mm, nn) = attenu_pat_rand(mm, nn) * exp(-1*s^2 / 6.5^2);
    end
end

signal_len = 2^13;         % length of signal, pad with 0
delta_k = delta_s / (2*pi);    % sampling period 
sampling_fre = 1 / delta_k;  % sampling frequency
delta_r = sampling_fre / signal_len;

perf_pat = zeros(signal_len, signal_len);  rand_pat = perf_pat;
perf_pat((signal_len/2+1-r_size) : (signal_len/2+1+r_size),...
              (signal_len/2+1-r_size) : (signal_len/2+1+r_size)) = attenu_pat_perf;
rand_pat((signal_len/2+1-r_size) : (signal_len/2+1+r_size),...
              (signal_len/2+1-r_size) : (signal_len/2+1+r_size)) = attenu_pat_rand;

pattern_len = 300;
pdf_perf = FourierHankel(perf_pat, signal_len, pattern_len);
pdf_rand = FourierHankel(rand_pat, signal_len, pattern_len);
pdf_perf = pdf_perf ./ max(pdf_perf,[],'all');
pdf_rand = pdf_rand ./ max(pdf_rand,[],'all');
rx = (-1*pattern_len : pattern_len) .* delta_r;
save('article_chapter2.1.mat');

%% figure 2.1 perfect align and random orientation
load('article_chapter2.1.mat');
fig = figure; set(gcf, 'position', [1200, 90, 500, 900]);    % subplot[5,2]. IT, IA, IM, sM and f(r) for each side line
ax(11) = axes('Position', [0.11, 0.83, 0.39, 0.145]);
ax(12) = axes('Position', [0.11, 0.642, 0.39, 0.145]);
ax(13) = axes('Position', [0.11, 0.454, 0.39, 0.145]);
ax(14) = axes('Position', [0.11, 0.265, 0.39, 0.145]);
ax(15) = axes('Position', [0.11, 0.05, 0.39, 0.145]);

ax(21) = axes('Position', [0.55, 0.83, 0.39, 0.145]);
ax(22) = axes('Position', [0.55, 0.642, 0.39, 0.145]);
ax(23) = axes('Position', [0.55, 0.454, 0.39, 0.145]);
ax(24) = axes('Position', [0.55, 0.265, 0.39, 0.145]);
ax(25) = axes('Position', [0.55, 0.05, 0.39, 0.145]);


% Perfect alignment 
imagesc(ax(11), sx, sy, IT_Pattern_perf);  colorbar(ax(11));  text(ax(11), -9, -9, '(A)', 'color', 'w');
ylabel(ax(11), 'sy ($\AA^{-1}$)','Interpreter','latex'), % xlabel(ax(11), 'sx ($\AA^{-1}$)','Interpreter','latex'),
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(11), 'FontSize', 10);  title(ax(11), 'Perfect Alignment-IT Pattern');

imagesc(ax(12), sx, sy, IA_Pattern_perf);  colorbar(ax(12));  text(ax(12), -9, -9, '(B)', 'color', 'w');
ylabel(ax(12), 'sy ($\AA^{-1}$)','Interpreter','latex'), %xlabel(ax(12), 'sx ($\AA^{-1}$)','Interpreter','latex'), 
% set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(12), 'FontSize', 10); title(ax(12), 'IA Pattern');

imagesc(ax(13), sx, sy, IM_Pattern_perf);  colorbar(ax(13));  text(ax(13), -9, -9, '(C)', 'color', 'w');
ylabel(ax(13), 'sx ($\AA^{-1}$)','Interpreter','latex'), %xlabel(ax(13), 'sy ($\AA^{-1}$)','Interpreter','latex')
% set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(13), 'FontSize', 10);  title(ax(13), 'Im Pattern');

imagesc(ax(14), sx, sy, M_Pattern_perf);  colorbar(ax(14));  text(ax(14), -9, -9, '(D)', 'color', 'w');
ylabel(ax(14), 'sx ($\AA^{-1}$)','Interpreter','latex'), xlabel(ax(14), 'sy ($\AA^{-1}$)','Interpreter','latex')
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(14), 'FontSize', 10); title(ax(14), 'sM Pattern'); 

pcolor(ax(15), rx, rx, pdf_perf);  colorbar(ax(15));  shading(ax(15), 'interp'), text(ax(15), -3.9, 3.9, '(E)', 'color', 'w');
xlabel(ax(15), 'rx ($\AA$)','Interpreter','latex'), ylabel(ax(15), 'ry ($\AA$)','Interpreter','latex')
set(gca, 'xtick', -3 : 3);  set(gca, 'xticklabel', [-3, -1, 0, 1, 3]); 
set(gca, 'ytick', -3 : 3);  set(gca, 'yticklabel', [-3, 1, 0, 1, 3]); 
set(ax(15), 'FontSize', 10);  title(ax(15), 'f(r) Pattern');

% random alignment
imagesc(ax(21), sx, sy, IT_Pattern_rand);  colorbar(ax(21));  text(ax(21), -9, -9, '(a)', 'color', 'w');
% xlabel(ax(21), 'sx ($\AA^{-1}$)','Interpreter','latex'), ylabel(ax(21), 'sy ($\AA^{-1}$)','Interpreter','latex')
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(21), 'FontSize', 10);  title(ax(21), 'Random Orientation-IT Pattern');

imagesc(ax(22), sx, sy, IA_Pattern_rand);  colorbar(ax(22));  text(ax(22), -9, -9, '(b)', 'color', 'w');
% xlabel(ax(22), 'sx ($\AA^{-1}$)','Interpreter','latex'), ylabel(ax(22), 'sy ($\AA^{-1}$)','Interpreter','latex')
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(22), 'FontSize', 10); title(ax(22), 'IA Pattern');

imagesc(ax(23), sx, sy, IM_Pattern_rand);  colorbar(ax(23));  text(ax(23), -9, -9, '(c)', 'color', 'w');
% xlabel(ax(23), 'sx ($\AA^{-1}$)','Interpreter','latex'), ylabel(ax(23), 'sy ($\AA^{-1}$)','Interpreter','latex')
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(23), 'FontSize', 10);  title(ax(23), 'Im Pattern');

imagesc(ax(24), sx, sy, M_Pattern_rand);  colorbar(ax(24));  text(ax(24), -9, -9, '(d)', 'color', 'w');
xlabel(ax(24), 'sx ($\AA^{-1}$)','Interpreter','latex'),% ylabel(ax(24), 'sy ($\AA^{-1}$)','Interpreter','latex')
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(ax(24), 'FontSize', 10); title(ax(24), 'sM Pattern'); 

pcolor(ax(25), rx, rx, pdf_rand);  colorbar(ax(25));shading(ax(25), 'interp'); colormap jet; text(ax(25), -3.9, 3.9, '(e)', 'color', 'w');
xlabel(ax(25), 'rx ($\AA$)','Interpreter','latex'),% ylabel(ax(25), 'ry ($\AA$)','Interpreter','latex')
set(gca, 'xtick', -3 : 3);  set(gca, 'xticklabel', [-3, -1, 0, 1, 3]); 
set(gca, 'ytick', -3 : 3);  set(gca, 'yticklabel', [-3, 1, 0, 1, 3]); 
set(ax(25), 'FontSize', 10);  title(ax(25), 'f(r) Pattern');