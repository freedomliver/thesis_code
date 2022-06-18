% Program for calculation of molecular alignment from Markus Guehr, 09Feb 2007
% Based on article Ortigoso et al, JCP 1999   Modified by Jinc 2022/3/31,
% Fixed at 2022/4/6
% 计算 N2, O2, I2, CS2, CO2分子的准直度(cos^2)时间演化(CF3I缺少极化参数)
% clc;  clear; close all;  
%{
epsilon0 = 8.854e-12;
e = 1.6022e-19;               % Elementary charge       [C]
k = 1.38066e-23;             % Boltzmann constant     [J/K]
c = 2.9979e8;                   % Speed of light in m/s   [L/S]
hbar = 1.054589e-34;      %                                       [J*S]

%% 常数和变量设置，更改分子信息（Mole，B0，delta_alpha，Boltzmann_statistics_ratio）
% 温度TemperatureK，以及激光参数（pulse_FWHM，Power，reprate，Laser_focus）
Mole = 'CF3I';
% Rotational constant in ground state in wavenumbers (1/cm)
% B0 = 1.998, 1.445, 0.037, 0.1091, 0.39, 1.52325;  for N2, O2, I2, CS2, CO2，CF3I
% B0=0.0508;% for CF3I
% B0=0.247;% for CH3I
B0 = 0.0506;                         %                                       [1/cm]
B = B0 * hbar * 2*pi * c * 1e2;    % Rotational constant in ground state in Joules   [J]

% Anisotropy of the polarizability in m^3 - but cgs units
% Delta_alpha = 0.696e-30, 1.1e-30, 6.69e-30, 9.6e-30, 2.1e-30;  cv
% For N2,O2, I2, CS2, CO2, CF3I
% delta_alpha=2e-30;% estimated for CF3I
% delta_alpha=2.5e-30;% estimated for CH3I
delta_alpha = 2e-30;                % For CO2   电离阈值约10^13 W/cm2   [m^3]

% For TemperatureK = 5 : 25 : 105,  for I = 0.1 : 0.2 : 1.1
TemperatureK = 10;                    % In Kelvin            [K]
Temperature = k * TemperatureK / B;             %      [1]

pulse_FWHM = 100e-15;            % FWHM duration of the intesity of the pulse  [S]
Power = 0.1;               % Power in watt                   [P]
reprate = 50;                % Rep rate in Hz                   [Hz]
pulse_energy = Power / reprate;% Pulse energy     [J]
disp(['Pulse energy is ' num2str(pulse_energy*1e6) 'uJ']);
Laser_focus = [100 100] * 3;        % Laser focus(FWHM)    [um]
Int_fac = 1;                   % Average intensity/Peak intensity
focus_area = Laser_focus(1)/2 * Laser_focus(2)/2 * pi / log(2)/Int_fac;   % log(2) from FWHM to 1/e. [um^2]

I = pulse_energy/pulse_FWHM/focus_area*1e12;  %  Convert into SI unit.      [W/m^2]
E_density =  pulse_energy/focus_area*1e8;            % [J/cm^2]
disp(['The Engry density used is ', num2str(E_density), ' J/cm^2']);
E0 = sqrt(I * 4*pi/c);    % Time averaged Electric field in cgs unit

%%
% Boltzmann statistics
addpath(genpath('D:\Program Files\MATLAB\cal\Alignment'));
global Jmax;
Jmax = 100;
pop = Boltzmann(Temperature);                 % Parameter is temperature given in units of B 
% Strength of the pulse in units of the rotational constant ***this is delta omega in the ortigoso paper
delta_omega = delta_alpha * E0^2 / B/2;   % [m^2]
EstimateJmax = double(int16(pulse_FWHM / (1/(B0*1e2)/c / delta_omega / 2/pi))); % How to get the Jmax

% Calculate highest initial J and population threshold
[qq] = max(pop);
disp(['Max population is ', num2str(qq), ', while the last population is ', num2str(pop(Jmax+1))]);
thre = qq*0.01;
ind = 1 : length(pop);
[qq2, Jhighpop] = max(ind(pop>thre));
Jmax = ceil((Jhighpop+EstimateJmax)/10) * 10;  % Set Jmax 10~20 higher than Estimated Jmax
upperJ = Jmax - 4;               % The highest rem. populated state due to temp
occ = zeros(upperJ, 2*upperJ+1, Jmax);  phase = occ;

disp(['Highest initial J occupied is  ', num2str(Jhighpop-1), '.']);
disp(['Estimate Highest excited J occupied is  ', num2str(Jhighpop-1+EstimateJmax), '.']);
disp(['Use Jmax =  ', num2str(Jmax), '.']);
disp(['After adjusted, the last population is ', num2str(pop(upperJ+1)), ', the ratio is ', num2str(pop(upperJ+1)./qq)]);

%%
% Propagation while pulse is on
% Sigma for the electric field gaussian as used in subroutine coefficients
global sigma; tic;
sigma1 = sqrt(pulse_FWHM^2 / (4*log(2)));    % [S]
sigma = sigma1 * B / hbar;         % Conver unit into hbar/B, consistent with paper  [1]
for Jstart = 0 : 1 : upperJ     % Jstart is the initially populated J - due to temperature
    Mdummy = 0;
    if mod(upperJ-Jstart, 10) == 0
        disp(['Remaining J # is  ', num2str(upperJ-Jstart+1), '.']);
        toc;
    end
    for Mstart = -Jstart : 1 : Jstart
        Mdummy = Mdummy + 1;
        % First parameter is parameter Delta omega, second is the initial J state populated, third the initial M
        c = coefficients(delta_omega, Jstart, Mstart);
        occ(Jstart+1, Mdummy, :) = abs(c);                 % Occupancy of state, its square normalized!!
        phase(Jstart+1, Mdummy, :) = angle(c);         % Phase of state
    end
end
toc;

%%
% Calculate the cos^2 expectation value  
t = 0 : 0.001 : 6.5;     % In unit of hbar/B
Jdummy = 0;  excitedJhighocc = 0; ciniJ = 0; ciniM = 0; 
for indexJtemp = 1 : upperJ   %this is the Jtemp+1
    Jtemp = indexJtemp - 1;
    population = pop(indexJtemp);
    if population >= thre
        Mdummy = 0;
        for indexMtemp = 1 : (2*(Jtemp) + 1)   %this is an index for M
            M = indexMtemp - 1 - Jtemp;
            % Part for abs calculation
            dummy = 0;
            for index = 1 : (Jmax-4)
                J = index - 1;
                if abs(occ(indexJtemp, indexMtemp, index)) > 0.05 && J > excitedJhighocc
                    excitedJhighocc = J;   ciniJ = Jtemp;   ciniM = M;
                end
                dummy = abs(occ(indexJtemp, indexMtemp, index)) .^ 2*cos2JM(J, M) + dummy;
            end
            timeindep = dummy;             % Diagonal terms in sum over J&J'
            % Part for phase calculation
            dummy = 0;
            for index = 1 : (Jmax-4)
                J = index - 1;
                phasep = phase(indexJtemp, indexMtemp, index) - phase(indexJtemp, indexMtemp, index+2);
                dummy = abs(occ(indexJtemp,indexMtemp,index)) * abs(occ(indexJtemp,indexMtemp,index+2))...
                    * (cos2Jp2M(J, M)) * cos((4*J+6).*t + phasep) + dummy;
                % 4J+6 is the energy difference between J+2 and J states in unit of B
            end
            timedep = dummy;                % Off diagonal terms in sum over J&J', only real part is considered

            Mdummy = 1/(2*Jtemp+1) * (2*timedep+timeindep) + Mdummy;  % Devide by M state number
        end
        Jdummy = population * Mdummy + Jdummy;
    end
end
toc;
time = t * hbar / B * 1e12;  cossquare = Jdummy;
figure;  plot(time, cossquare, '-x');  grid on;  xlabel('time(ps)');  ylabel('<cos^2\theta>');
title([Mole, ' with E-density = ', num2str(E_density), 'J/cm^2 at Temperature = ', num2str(TemperatureK), 'K']);
% disp(['Highest excited J is ', num2str(excitedJhighocc),', corresponding initial J is ' ,num2str(ciniJ), ...
%     ', corresponding initial M is ',num2str(ciniM),'.']);        % Largest J excited
% disp(['Laser intensity peaked at ', num2str(3*1e15*sigma1),' fs before t=0.']); 

% filenameI = num2str(E_density);
% filenameT = num2str(TemperatureK);
% filename = strcat('Cossquare-', Mole,'-Engery_density_', filenameI, '-Temp_', filenameT, 'K.dat');
% savevar(:, 1) = time;
% savevar(:, 2) = cossquare;
% save(filename, 'savevar', '-ASCII');

%%
[m1, p1] = max(cossquare);
[m2, p2] = min(cossquare);
t2 = 0.001 .* (0 + p1);                                         %选取计算最大还是最小
resolution = [1, 41];  wavepacketsq = 0;  wavepacket = 0;
for indexJtemp = 1 : upperJ                             %this is the Jtemp+1
    for indexMtemp = 1 : (2*(indexJtemp-1) + 1)  %this is an index for M
        Jtemp = indexJtemp - 1;
        M = indexMtemp - 1 - Jtemp;
        population = pop(indexJtemp) / (2*Jtemp+1);  
       %now begins the coherent superposition of wavefunctions for each initially populated J level
        wavepacket = 0;
        for index = 1 : Jmax
            J = index - 1;
            if occ(indexJtemp, indexMtemp, index) ~= 0
                wavepacket = occ(indexJtemp, indexMtemp, index) * spharm(J, M, [resolution(1), resolution(2)], 1)...
                                        * exp(-1i*J*(J+1)*t2 + 1i*phase(indexJtemp, indexMtemp, index)) + wavepacket;
            end
        end
        wavepacketsq = population * abs(wavepacket).^2 + wavepacketsq;
    end
end

save('2022_5_17.mat');
THETA = linspace(0, pi, resolution(2)); distrib = sum(wavepacketsq, 2);
% THETA = THETA(1 : 41); distrib = distrib(1 : 41);
% THETA = reshape(THETA, length(THETA), 1); distrib = reshape(distrib, length(distrib), 1);
% distrib = distrib ./ sum(distrib' .* sin(THETA)); 
% 为什么要除以sin的归一化？是为了方便计算衍射班？
figure;
plot(THETA, distrib, '-x'); grid on;
xlabel('theta(θ)'); ylabel('归一化分布');
title(['Maximum cos^2 of ', Mole, ' is ',num2str(m1), ', Temperature is ', num2str(TemperatureK), 'K',...
        newline, 'The maximun time is ', num2str(t2*hbar/B*1e12), 'ps after laser aligned']);
% savedistrib(:, 1) = THETA;
% savedistrib(:, 2) = distrib;
% filename = strcat('Max-distrib-', Mole, '-Engery_density_', filenameI, '-Temp_', filenameT, 'K.dat');
% save(filename, 'savedistrib', '-ASCII');
%}
%% figure 2.2 Alignment evolution of  CF3I over time
load('CF3I_Alignment_simulation.mat');
figure; set(gcf, 'position', [900, 90, 900, 400]);    
ax(11) = axes('Position', [0.08, 0.15, 0.37, 0.74]);
ax(12) = axes('Position', [0.57, 0.15, 0.38, 0.74]);
% ax(2) = axes('Position', [0.1, 0.12, 0.85, 0.33]);

plot(ax(11), time(1:5000), cossquare(1:5000), '-');  grid(ax(11), 'on');  
xlabel(ax(11), 'time(ps)');  ylabel(ax(11), '<cos^2\theta>');  text(ax(11), 50, 0.55, '(A)', 'color', 'k');
title(ax(11), [Mole, ' with ', num2str(E_density), 'J/cm^2 at Temperature = ', num2str(TemperatureK), 'K']);
set(ax(11), 'FontSize', 14); 

plot(ax(12), THETA, distrib, '-x'); grid(ax(12), 'on');
% plot(ax(12), THETA, Evolution(p2, :), '-x'); grid(ax(12), 'on');
xlabel(ax(12), 'theta(θ)'); ylabel(ax(12), '归一化分布');  text(ax(12), 0.1, 0.14, '(B)', 'color', 'k');
title(ax(12), ['Max-cos^2  = ', num2str(m1), ', T = ', num2str(TemperatureK), 'K', ]);
set(ax(12), 'FontSize', 14); 

% imagesc(ax(2), time, THETA, Evolution');    text(ax(2), 0.4, 0.2, '(C)', 'color', 'k');
% xlabel(ax(2), 'Time delay (ps)'); ylabel(ax(2), 'The entire angle'); colorbar(ax(2));
% set(gca, 'FontSize', 14);  


%{
%%  乘子算法的效果
addpath(genpath('D:\Program Files\MATLAB\cal\src'));    % Add the path of program package
DPWAfile = 'D:\Program Files\MATLAB\cal\src\DPWA\';   % Diffraction amplitude path (DPWA file)
sx = 0 : 0.0264 : 10;  sy = -10 : 0.0264 : 10;
smap = zeros(length(sy), length(sx));      % scattering amplitude
for ii = 1 : length(sy), for jj = 1 : length(sx), smap(ii, jj) = norm([sy(ii), sx(jj)]); end; end
[fC, ~] = DPWAcpu(smap, 6, 3000, DPWAfile);
[fF, ~ ] = DPWAcpu(smap, 9, 3000, DPWAfile);
[fI, ~ ] = DPWAcpu(smap, 53,3000, DPWAfile);

rC = [0; 0; 0];  rI = [0; 2.1655; 0];    % C at the orign, I on the axis y
rF1 = [1.2405; -0.4630; 0];
rF2 = [-0.6203; -0.4630; 1.0743];
rF3 = [-0.6203; -0.4630; -1.0743];
molecular = [rC, rI, rF1, rF2, rF3];

cossquare = load('D:\Program Files\MATLAB\cal\thesis\Cossquare-CF3I-Engery_density_1.9612-Temp_10K.dat');
Evolution = load('D:\Program Files\MATLAB\cal\thesis\Evolution-CF3I-Engery_density_1.9612-Temp_10K.dat');
[m1, p1] = max(cossquare);  [m2, p2] = min(cossquare);
distri_max = Evolution(p1(2), :); distri_min = Evolution(p2(2), :);
beta = linspace(0, pi, length(distri_max));

alpha = 0 : (2*pi/40) : 2*pi;   gamma = 0 : (2*pi/40) : 2*pi;  
AngularDis_min = distri_min .* sin(beta);       
AngularDis_max = distri_max .* sin(beta);       
AngularDis_rand = sin(beta);
AngularDis_hmax = Evolution(1528, :) .* sin(beta);

% beta = 0;  AngularDis = 1;
% 
% figure;  plot(cossquare(:,2), '-x');  grid on;  xlabel('time(ps)');  ylabel('<cos^2\theta>');
% title([Mole, ' with E-density = ', num2str(E_density), 'J/cm^2 at Temperature = ', num2str(TemperatureK), 'K']);

% figure, plot(beta, AngularDis_max, '-x'); grid on
% xlabel('$\theta$', 'Interpreter', 'latex'),  ylabel('f($\theta$) (arb.unit)', 'Interpreter', 'latex');
% set(gca, 'linewidth', 1),  set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
% set(gca, 'Xtick', [0,pi/2,pi]),  set(gca, 'XTicklabel', {'0','\pi/2','\pi'});
%%
AngularDis = AngularDis_hmax;
IT = smap .* 0;  IA = smap .* 0;  tic
for beta_index = 1 : length(beta) 
    theta = beta(beta_index);
    layer1 = repmat(sx, length(sy), 1);       % sx
    layer2 = repmat(sy', 1, length(sx));      % sy
    layer3 = zeros(length(sy), length(sx)); % sz=0
    s = cat(3, layer1, layer2, layer3);          % dimension3: [sx,sy,sz]
    
    It = smap .* 0;  Ia = smap .* 0;
    for alpha_index = 1 : length(alpha)
        phi = alpha(alpha_index);       
        for gamma_index = 1 : length(gamma) % gamma
            roll = gamma(gamma_index);
            r = RotationMatrixYZY(phi, theta, roll) * molecular; % r
            diffraciton = fC .* exp(1i.*sum(s.*reshape(r(:,1), [1,1,3]), 3)) + ...
                                  fI .* exp(1i.*sum(s.*reshape(r(:,2), [1,1,3]), 3)) + ...
                                  fF .* exp(1i.*sum(s.*reshape(r(:,3), [1,1,3]), 3)) + ...
                                  fF .* exp(1i.*sum(s.*reshape(r(:,4), [1,1,3]), 3)) + ...
                                  fF .* exp(1i.*sum(s.*reshape(r(:,5), [1,1,3]), 3));                % sum(f*exp(isr))
            It = It + abs(diffraciton).^2;         
            Ia = Ia + abs(fC).^2 + abs(fI).^2 + abs(fF).^2.*3;          
        end
    end
    DiffractionTheta = It ./ (length(alpha) * length(gamma));                          % single molecule
    Ia_theta = Ia ./ (length(alpha) * length(gamma));
    
    if theta - pi < 0.00001  % theta=pi
        IT = IT + DiffractionTheta .* AngularDis(beta_index) ./ 2; % for symmetry
        IA = IA + Ia_theta .* AngularDis(beta_index) ./ 2;
    else
        IT = IT + DiffractionTheta .* AngularDis(beta_index);
        IA = IA + Ia_theta .* AngularDis(beta_index);
    end 
    toc
end
% I total
IA = IA ./ sum(AngularDis);  IT = IT ./ sum(AngularDis);

IT_Pattern = [fliplr(IT(:, 1:end)), IT] ./ max(IT, [], 'all');
IA_Pattern = [fliplr(IA(:, 1:end)), IA] ./ max(IA, [], 'all');
IM = IT - IA;  IM_Pattern = [fliplr(IM(:, 1:end)), IM] ./ max(IT, [], 'all');
M = IM ./ IA;  M_left = fliplr(M(:, 1:end));  M_Pattern = [M_left, M];
IT_Pattern_hmax = IT_Pattern;
IA_Pattern_hmax = IA_Pattern;
% % combine = 最大准直度和只有最大准直度一半位置的差值 + 模拟的随机IM
% IM_combine = M_Pattern_max_hmax +IM_Pattern_rand;  
% % 用IM_combine 和模拟的随机IA 计算约化衍射图sM
% M_combine = IM_combine ./ IA_Pattern;
% % 用M_combine 通过傅里叶汉克而变换得到PDF
figure
pcolor(sy,sy,M_Pattern), colorbar; shading interp,colormap jet
xlabel('sx ($\AA^{-1}$)','Interpreter','latex'), ylabel('sy ($\AA^{-1}$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',20)
set(gca,'XTick',-12:4:12),set(gca,'XTicklabel',-12:4:12)
set(gca,'YTick',-12:4:12),set(gca,'YTicklabel',-12:4:12)

% M_Pattern_rand = M_Pattern;
% IM_Pattern_rand = IM_Pattern;
%}
%%  画论文用的图
load('chengzisuanfa.mat');
figure; set(gcf, 'position', [400, 300, 1000, 300]);
ax(1) = axes('Position', [0.05, 0.19, 0.27, 0.72]);

pcolor(ax(1),rx,rx,pdf_max), colorbar(ax(1)); shading(ax(1), 'interp'), colormap jet
xlabel('rx ($\AA$)','Interpreter','latex'), ylabel('ry ($\AA$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',13)
set(gca,'xtick',-3:3), set(gca,'ytick',-3:3)
title('PDF-maximun cos^2 \theta');

ax(2) = axes('Position', [0.37, 0.19, 0.27, 0.72]);
pcolor(ax(2),rx,rx,pdf_hmax), colorbar(ax(2)); shading(ax(2), 'interp'), colormap jet
xlabel('rx ($\AA$)','Interpreter','latex'), ylabel('ry ($\AA$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',13)
set(gca,'xtick',-3:3), set(gca,'ytick',-3:3)
title('PDF-half-max cos^2 \theta');

ax(3) = axes('Position', [0.69, 0.19, 0.27, 0.72]);
pcolor(ax(3),rx,rx,pdf_max_hmax), colorbar(ax(3)); shading(ax(3), 'interp'), colormap jet
xlabel('rx ($\AA$)','Interpreter','latex'), ylabel('ry ($\AA$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',13)
set(gca,'xtick',-3:3), set(gca,'ytick',-3:3)
title('PDF-of   \Delta sM');

%%
figure; set(gcf, 'position', [400, 300, 1000, 300]);
ax(1) = axes('Position', [0.05, 0.19, 0.27, 0.72]);

pcolor(ax(1),rx,rx,pdf_pri_1), colorbar(ax(1)); shading(ax(1), 'interp'), colormap jet
xlabel('rx ($\AA$)','Interpreter','latex'), ylabel('ry ($\AA$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',13)
set(gca,'xtick',-3:3), set(gca,'ytick',-3:3); text(ax(1), -3.9, 3.9, '(A)', 'color', 'w');
title('PDF with ef = 1');

ax(2) = axes('Position', [0.37, 0.19, 0.27, 0.72]);
pcolor(ax(2),rx,rx,pdf_pri_05), colorbar(ax(2)); shading(ax(2), 'interp'), colormap jet
xlabel('rx ($\AA$)','Interpreter','latex'), ylabel('ry ($\AA$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',13)
set(gca,'xtick',-3:3), set(gca,'ytick',-3:3); text(ax(2), -3.9, 3.9, '(B)', 'color', 'w');
title('PDF with ef = 0.5');

ax(3) = axes('Position', [0.69, 0.19, 0.27, 0.72]);
pcolor(ax(3),rx,rx,pdf_pri_02), colorbar(ax(3)); shading(ax(3), 'interp'), colormap jet
xlabel('rx ($\AA$)','Interpreter','latex'), ylabel('ry ($\AA$)','Interpreter','latex')
set(gca,'linewidth',1), set(gca,'FontName','Times New Roman','FontSize',13)
set(gca,'xtick',-3:3), set(gca,'ytick',-3:3);   text(ax(3), -3.9, 3.9, '(C)', 'color', 'w');
title('PDF with ef = 0.2');


%%    处理准直实验和模拟结果  主要要改变模拟计算的文件夹
% This script is used for CF3I aligning experiment, the data used is 2021-8/4 and 2021-8/5
clc; clear;
addpath(genpath('D:\Program Files\MATLAB\cal\src'));   % Add the path for functions used

%%  Reading the experiment data
[PatternNames, DelayTime, numFolder, numPattern] = LoadingFolder('D:\气态');
TimeDelay = mean(DelayTime);                                           % calculating time delay
TimeDelay=2*(TimeDelay-TimeDelay(1)) ./ 1e7 ./ 3e8 .* 1e15;

BackFile = 'D:\气态\2021-08-04 CF3I准直\back.spe';       % loading background  
BackPatterns = ReadSPE(BackFile);

%% find hole center  寻找衍射图的中心位置，中心圆圈只和荧光屏有关系，不会随着电子抖动发生改变。
SinglePattern = double(imread(PatternNames{1,1}));
[numRowPix, numColPix] = size(SinglePattern);   
SizePix = [numRowPix, numColPix];                                     % 用一个数据获得衍射图大小

HolePattern = zeros(numRowPix, numColPix, numFolder);% 每个文件夹第一张图
for ii = 1 : numFolder,  HolePattern(:, :, ii) = double(imread(PatternNames{ii, 1}));  end 
HoleFinderPattern = median(HolePattern, 3);                      % Averaging all patterns，第一张图取平均值
HoleStrengthMax = 550;  HoleDeviation = 130;   
[HoleCenterRow, HoleCenterCol, HoleRadius] = FindHoleCenter(HoleFinderPattern, HoleStrengthMax, HoleDeviation);
HoleCenter = [HoleCenterRow, HoleCenterCol];                  % 洞中心的坐标

%% creating a single mask pattern and the datas in the hole are 0
MaskPatternSingle = ones(numRowPix, numColPix);          % 中间圆区域归为零
for ii = 1 : numRowPix
for jj = 1 : numColPix
if sqrt((ii-HoleCenterRow)^2+(jj-HoleCenterCol)^2) <= (HoleRadius+1), MaskPatternSingle(ii, jj) = 0; end; end
end
DarkCurrentFig = figure;               
pcolor(median(HolePattern, 3)),  shading interp,  colormap(jet),  caxis([480, 3000]) 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
set(gcf, 'position', [100, 100, 800, 800]);  axis equal

Flag = 0;         
Choice = questdlg('Do you want to remove dark current', 'Removing Dark Current', 'Yes', 'No', 'Yes');
if strcmp(Choice, 'Yes') == 1,  Flag = 1;  end
while Flag == 1
    DarkCurrent = getrect(DarkCurrentFig);
    DarkColumn = round(DarkCurrent(1)) : round(DarkCurrent(1))+round(DarkCurrent(3));
    DarkRow = round(DarkCurrent(2)) : round(DarkCurrent(2))+round(DarkCurrent(4));
    MaskPatternSingle(DarkRow, DarkColumn) = 0;
    Choice = questdlg('Do you want to continue', 'Removing Dark Current', 'Yes', 'No', 'No');
    if strcmp(Choice, 'No') == 1,  Flag = 0;  end     
end
%% circle fitting test
FindCenterRadius = 180;  RadiusCoefficient = [0.6, 1.6]; 
[MeanRegion, FitRegion] = PreFindPatternCenter(FindCenterRadius, RadiusCoefficient, SizePix, HoleCenter);

CircleThreshold = [0.95, 1.2];                                               % CircleThreshold: data range for fitting
[CenterRow,CenterCol,CenterRadius,TestPattern] =...
    FindPatternCenter(MeanRegion,FitRegion,medfilt2(SinglePattern,[9,9]),CircleThreshold,MaskPatternSingle);

figure;  % show circle fitting
imshow(TestPattern)
hold on
viscircles([HoleCenterCol, HoleCenterRow], HoleRadius, 'EdgeColor', 'r');
viscircles([HoleCenterCol, HoleCenterRow], FindCenterRadius*RadiusCoefficient(1), 'EdgeColor', 'g');
viscircles([HoleCenterCol, HoleCenterRow], FindCenterRadius*RadiusCoefficient(2), 'EdgeColor', 'g');
viscircles([CenterCol, CenterRow], CenterRadius, 'EdgeColor', 'b');

%% acquiring average patterns with different time delay
numPix = 951;  Range = (numPix-1)/2; 
TimePattern = zeros(numPix, numPix, numPattern); TimeMask = TimePattern;  TimeStd = TimePattern;
for ii = 1 : numPattern                                                         % processing the patterns with different time delay
    disp(['===being processed with different time delay: ',num2str(round(ii/numPattern*100)),'%===']);
    Pattern = zeros(numRowPix, numColPix, numFolder);  % initial variate to store patterns with same time delay
    MaskPattern = repmat(MaskPatternSingle, 1, 1, numFolder);
    ModifiedPattern = zeros(numPix, numPix, numFolder);% initial variate to store modified patterns 
    ModifiedMask = zeros(numPix, numPix, numFolder);   % with same time delay
    
    RegionRemoved = 9;                                                        % the size of removed region
    for jj = 1 : numFolder                                                        % processing the patterns with same time delay
        Pattern(:, :, jj) = double(imread(PatternNames{jj,ii}));%Pattern(:, :, jj) = ReadSIF_single(PatternNames{jj, ii});
                                                                                               % hot pixel removal(caused by X-Ray at EMCCD)
        [MaskPattern(:, :, jj), ~] = HotPixelRemoval(Pattern(:, :, jj), MaskPattern(:, :, jj), RegionRemoved);
        
        [CenterRow, CenterCol, ~, ~] =...                                 % find the center of the diffraciton pattern
            FindPatternCenter(MeanRegion, FitRegion, medfilt2(Pattern(:,:,jj),[9,9]), CircleThreshold, MaskPattern(:, :, jj));
%         CenterRow = 522;   % warning: temporary application
%         CenterCol = 443;   可以手动标中心，程序不好的话

        % align the patterns and cut the patterns to suitable size
        RowRange = (CenterRow-Range) : (CenterRow+Range);
        ColRange = (CenterCol-Range) : (CenterCol+Range);
        ModifiedPattern(:, :, jj) = Pattern(RowRange, ColRange, jj);
        ModifiedMask(:, :, jj) = MaskPattern(RowRange, ColRange, jj);
    end
    % calculating the mean pattern of the same delay and removing the hot pixels
    numCycles = 3;   ThresholdStd = 3;
    [TimePattern(:,:,ii), TimeMask(:,:,ii), TimeStd(:,:,ii)] = ...
                    AveragePattern(ModifiedPattern, ModifiedMask, numCycles, ThresholdStd);

    % Baseline subtraction (in the conner of the pattern without phosphor screen)
    BaselineWidth = 60;
    Baseline = sum(sum(TimePattern(end-BaselineWidth : end, end-BaselineWidth : end, ii)...
                                    .* TimeMask(end-BaselineWidth : end, end-BaselineWidth : end, ii)));
    Baseline = Baseline / sum(sum(TimeMask(end-BaselineWidth:end, end-BaselineWidth:end, ii)));
    TimePattern(:,:,ii) = TimePattern(:,:,ii) - Baseline;
end
tmp = TimePattern;  
for ii = 1 : numPattern,  TimePattern(:, :, ii) = tmp(:, :, numPattern-ii+1);  end
disp('===Patterns with different time delay have been acquired===')

figure  % test the quality of the average
subplot(1,3,1)
surf(TimePattern(:,:,1)), shading interp, colormap(jet), caxis([-10,2000])
subplot(1,3,2)
pcolor(TimeMask(:,:,1)), shading interp, colormap(jet)
subplot(1,3,3)
pcolor(TimeStd(:,:,1)), shading interp, colormap(jet)
set(gcf,'unit','centimeters','position',[1,10,50,12])

% load('8_5_Scan13_read.mat');
BackFile = 'D:\气态\2021-08-04 CF3I准直\back.spe';       % loading background  
BackPatterns = ReadSPE(BackFile);

%% Averaging back pattern
numPix = 951;  Range = (numPix-1)/2; 
RowRange = (CenterRow-Range) : (CenterRow+Range);
ColRange = (CenterCol-Range) : (CenterCol+Range);
BackPatterns = BackPatterns(RowRange, ColRange, :);
PreBackMask = MaskPatternSingle(RowRange, ColRange);
back_num = size(BackPatterns);
PreBackMask = repmat(PreBackMask, 1, 1, back_num(3));

[Back, BackMask, BackStd] = AveragePattern(BackPatterns, PreBackMask, 3, 1.5);
Baseline = sum(sum(Back(end-BaselineWidth:end, end-BaselineWidth:end) ...
                        .* BackMask(end-BaselineWidth:end, end-BaselineWidth:end)));
Baseline = Baseline / sum(sum(BackMask(end-BaselineWidth:end, end-BaselineWidth:end)));
Background = Back - Baseline;
disp('===Background pattern has been acquired===')
figure
surf(Background), shading interp, colormap(jet), caxis([-20,300])


%% remove background
TimePattern = TimePattern - Background;
TimeMask = TimeMask .* BackMask;

% medfilt2
MedFiltRegion = [9, 9];
for ii = 1 : numPattern
    TimePattern(:,:,ii) = medfilt2(TimePattern(:,:,ii), MedFiltRegion);
    TimeMask(:,:,ii) = medfilt2(TimeMask(:,:,ii) ./ TimeMask(:,:,ii), MedFiltRegion); 
end
TimeMask(isnan(TimeMask)) = 0;

NormRegion = 200 : 700;  NormFactor = zeros(1, numPattern);   % Normlization
for ii = 1 : numPattern,  NormFactor(ii) = nansum(nansum(TimePattern(NormRegion,NormRegion,ii)));  end
NormFactor = NormFactor ./ NormFactor(1);
for ii = 1 : numPattern,  TimePattern(:,:,ii) = TimePattern(:,:,ii) ./ NormFactor(ii);  end

figure
surf(TimePattern(:,:,end)), shading interp, colormap(jet), caxis([-20,700])

%% 计算模拟各项异性值  clc;clear;
% 通过PSF计算衍射信号的演化
cossquare = load('D:\Program Files\MATLAB\cal\thesis_code\Cossquare-CF3I-Engery_density_1.9612-Temp_10K.dat');
Evolution = load('D:\Program Files\MATLAB\cal\thesis_code\Evolution-CF3I-Engery_density_1.9612-Temp_10K.dat');
figure; 
plot(cossquare, '-x');

addpath(genpath('D:\Program Files\MATLAB\cal\src'));    % Add the path of program package
DPWAfile = 'D:\Program Files\MATLAB\cal\src\DPWA\';   % Diffraction amplitude path (DPWA file)
sx = 0 : 0.025 : 10;  sy = -10 : 0.025 : 10;
smap = zeros(length(sy), length(sx));      % scattering amplitude
for ii = 1 : length(sy), for jj = 1 : length(sx), smap(ii, jj) = norm([sy(ii), sx(jj)]); end; end
[fC, ~] = DPWAcpu(smap, 6, 3000, DPWAfile);
[fF, ~ ] = DPWAcpu(smap, 9, 3000, DPWAfile);
[fI, ~ ] = DPWAcpu(smap, 53,3000, DPWAfile);

rC = [0; 0; 0];  rI = [0; 2.1655; 0];    % C at the orign, I on the axis y
rF1 = [1.2405; -0.4630; 0];
rF2 = [-0.6203; -0.4630; 1.0743];
rF3 = [-0.6203; -0.4630; -1.0743];
molecular = [rC, rI, rF1, rF2, rF3];

%计算完美准直的信号，需要重新计算，按照荧光屏的转移动量大小计算
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
        for gamma_index = 1 : length(gamma) % gamma
            roll = gamma(gamma_index);
            r = RotationMatrixYZY(phi, theta, roll) * molecular; % r
            diffraciton = fC .* exp(1i.*sum(s.*reshape(r(:,1), [1,1,3]), 3)) + ...
                                  fI .* exp(1i.*sum(s.*reshape(r(:,2), [1,1,3]), 3)) + ...
                                  fF .* exp(1i.*sum(s.*reshape(r(:,3), [1,1,3]), 3)) + ...
                                  fF .* exp(1i.*sum(s.*reshape(r(:,4), [1,1,3]), 3)) + ...
                                  fF .* exp(1i.*sum(s.*reshape(r(:,5), [1,1,3]), 3));                % sum(f*exp(isr))
            It = It + abs(diffraciton).^2;         
            Ia = Ia + abs(fC).^2 + abs(fI).^2 + abs(fF).^2.*3;          
        end
    end
    IT = IT + It ./ (length(alpha) .* length(gamma));
    IA = IA+ Ia ./ (length(alpha) .* length(gamma));
end
IT_Pattern = [fliplr(IT(:, 2:end)), IT] ./ max(IT, [], 'all');
IA_Pattern = [fliplr(IA(:, 2:end)), IA] ./ max(IA, [], 'all');
IM = IT - IA;  IM_Pattern = [fliplr(IM(:, 2:end)), IM] ./ max(IT, [], 'all');
M = IM ./ IA;  M_left = fliplr(M(:, 2:end));  M_Pattern = [M_left, M];

%% 测试 对于完美准直衍射信号，沿着几乎完美准直角分布PSF  通过PSF计算得到的衍射信号
time_some = 1:10:3550;
GET_CF3I_pattern_diff = zeros(pat_size, pat_size,length(time_some));  tic
% IM_PATTERN_CF3I_rand = getPSFpattern(IM_Pattern, ones(1,41), psf);
for i = 1 : length(time_some)
    angulaDis_t = time_some(i); 
    GET_CF3I_pattern_diff(:,:, i)=getPSFpattern(IM_Pattern, Evolution(angulaDis_t,:), psf);% - IM_PATTERN_N2_rand; 
    toc 
end

save('3.4_CF3I_anistropy_exper_simu.mat');

%%  模拟计算各项异性值的部分
cc = 201;   bias = 101;  a1 = 59;  b1 = 31;    
x1 = cc-bias-a1;  x2 = cc-bias+a1; x3 = cc+bias-a1;  x4 = cc+bias+a1;   y1 = cc-b1; y2 = cc+b1;

diff0 = zeros(1, length(time_some));
for i=1:length(time_some)
    figure(97);  imagesc(GET_CF3I_pattern_diff(:, :, i)); title(['i=',num2str(i)]); pause(0.001);   colorbar;  hold on 
    plot([ x1,x1,x2,x2,x3,x3,x4,x4,y1,y2,y1,y2,y1,y2,y1,y2],...
            [y1,y2,y1,y2,y1,y2,y1,y2,x1,x1,x2,x2,x3,x3,x4,x4], 'xr', 'MarkerSize', 15);
    
    ave_a = mean(GET_CF3I_pattern_diff(x1:x2, y1:y2, i), 'all');
    ave_b = mean(GET_CF3I_pattern_diff(x3:x4, y1:y2, i),'all');
    ave_c = mean(GET_CF3I_pattern_diff(y1:y2, x1:x2, i),'all');
    ave_d = mean(GET_CF3I_pattern_diff(y1:y2, x3:x4, i), 'all');
    diff0(i) = - (ave_a + ave_b - ave_c - ave_d);
end

figure;
plot(time_some, diff0);
%}
%%  CF3I各向异性模拟图
% load('CF3I_experiment_simulation_anisotropy.mat');
c = 2.9979e8;  hbar = 1.054589e-34;  B0=0.0508;  B = B0 * hbar * 2*pi * c * 1e2; 
time = time_some * hbar / B * 1e9;
align_t = cossquare(:,1);  align_cossquare = cossquare(:,2);
figure;
yyaxis left; plot(time, diff0); ylabel('anisotropy'); hold on; yyaxis right;
plot(timeo+3.44, (diff_exp ./ 100),'o-r');ylabel('<cos^2\theta>');
xlabel('times(ps)');title('氮气分子准直度和各向异性模拟图'); set(gca, 'FontSize', 13);
axes('position', [0.20, 0.18, 0.18,0.25]);
imagesc(GET_CF3I_pattern_diff(:, :, 300));  hold on
plot([ x1,x1,x2,x2,x3,x3,x4,x4,y1,y2,y1,y2,y1,y2,y1,y2],...
            [y1,y2,y1,y2,y1,y2,y1,y2,x1,x1,x2,x2,x3,x3,x4,x4], 'xr', 'MarkerSize', 4);
set(gca, 'xtick', [], 'xticklabel', []);set(gca, 'ytick', [], 'yticklabel', []);

%%  CF3I各向异性模拟图
c = 2.9979e8;  hbar = 1.054589e-34;  B0=0.0508;  B = B0 * hbar * 2*pi * c * 1e2; 
time = time_some * hbar / B * 1e9;
align_t = cossquare(:,1);  align_cossquare = cossquare(:,2);
figure;
yyaxis left; plot(time, diff0); ylabel('anisotropy-simulation'); hold on; yyaxis right;
plot(timeo+2.5, -(diff_exp ./ 10),'o-r');ylabel('anisotropy-experiment');
xlabel('times(ps)');title('CF3I分子各向异性模拟与实验图'); set(gca, 'FontSize', 13);

axes('position', [0.18, 0.18, 0.19,0.3]);
yyaxis left; plot(time(1:11), diff0(1:11)); set(gca, 'ytick', [], 'yticklabel', []);
hold on; yyaxis right; plot(timeo(1:24)+2.8, -(diff_exp(1:24) ./ 100),'x-r');
set(gca, 'xtick', [], 'xticklabel', []);set(gca, 'ytick', [], 'yticklabel', []);

axes('position', [0.45, 0.18, 0.26,0.3]);
yyaxis left; plot(time(307:324), diff0(307:324)); set(gca, 'ytick', [], 'yticklabel', []);
hold on; yyaxis right;  plot(timeo(25:44)+2.8, -(diff_exp(25:44) ./ 100),'x-r');
set(gca, 'xtick', [], 'xticklabel', []);set(gca, 'ytick', [], 'yticklabel', []);


% load('article_chapter3.4.2_CF3I_anistropy.mat');

%%  计算实验各项异性值的部分
cc = 201;   bias = 101;  a1 = 59;  b1 = 31;    
x1 = cc-bias-a1;  x2 = cc-bias+a1; x3 = cc+bias-a1;  x4 = cc+bias+a1;   y1 = cc-b1; y2 = cc+b1;
time_pattern_size = size(TimePattern);
% backsize = size(IM_PATTERN_CF3I_rand);
diff_exp = zeros(1, time_pattern_size(3));
GET_EXPCF3I_pattern = zeros(pat_size, pat_size, time_pattern_size(3));
for i=1:time_pattern_size(3)
    GET_EXPCF3I_pattern(:, :, i) = TimeDiff(((length(TimePattern)-400)/2+1) : (length(TimePattern)+400)/2, ...
    ((length(TimePattern)-400)/2+1) : (length(TimePattern)+400)/2, i);
    figure(97);  imagesc(GET_EXPCF3I_pattern(:, :, i)); title(['i=',num2str(i)]); pause(0.001);   colorbar;  hold on 
    plot([ x1,x1,x2,x2,x3,x3,x4,x4,y1,y2,y1,y2,y1,y2,y1,y2],...
            [y1,y2,y1,y2,y1,y2,y1,y2,x1,x1,x2,x2,x3,x3,x4,x4], 'xr', 'MarkerSize', 15);
    
    ave_a = mean(GET_CF3I_pattern_diff(x1:x2, y1:y2, i), 'all');
    ave_b = mean(GET_CF3I_pattern_diff(x3:x4, y1:y2, i),'all');
    ave_c = mean(GET_CF3I_pattern_diff(y1:y2, x1:x2, i),'all');
    ave_d = mean(GET_CF3I_pattern_diff(y1:y2, x3:x4, i), 'all');
    diff_exp(i) = - (ave_a + ave_b - ave_c - ave_d);% ./ (ave_a + ave_b + ave_c + ave_d) ;
end

figure;
plot(time_some(1:44), diff_exp);
diff_exp = diff0;
% diff = fliplr(diff0) - mean(diff0(end-6+1:end-1+1));
[stm,stt]=max(diff0(1:end));
% tt=fliplr(t);
% time0=tt-tt(stt); %把峰值设为0点
timeo = (TimeDelay - stt.*(TimeDelay(2) - TimeDelay(1))) ./ 1000;
figure;plot(timeo, (diff_exp ./ 100),'o-r');
% hold on; plot(time(1:3200)-2.613, cossquare(1:3200) , '-x');  grid on; 
% axis([-3 333 -1 1]);
breakxaxis([20  321]);
xlabel('t (ps)'); 
ylabel('anisotropy (%)');  
set(gca,'FontSize',20);
% figure;

