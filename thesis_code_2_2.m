% This part is prepared for undergraduate thesis chapter 2.2
% This script is taking about N2 alignment
clc;  clear;  
addpath(genpath('D:\Program Files\MATLAB\cal\Alignment'));  
Mole = 'N2';
Mole1 = 'I2';
Mole2 = 'CS2';
% Some constents
epsilon0 = 8.854e-12;  e = 1.6022e-19;            
k = 1.38066e-23;  c = 2.9979e8;  hbar = 1.054589e-34;    

% 这里是CS2
B0 = 0.1091;  B = B0 * hbar * 2*pi * c * 1e2;  delta_alpha = 9.6e-30; 
TemperatureK = 30;  Temperature = k * TemperatureK / B;      
pulse_FWHM = 100e-15; Power = 0.08; reprate = 50;           
pulse_energy = Power / reprate;  focus_area = 150^2 * pi / log(2);   
I = pulse_energy/pulse_FWHM/focus_area*1e12;  E0 = sqrt(I * 4*pi/c);  
E_density =  pulse_energy/focus_area*1e8;   
disp(['The Engry density used is ', num2str(E_density), ' J/cm^2']);

%
global Jmax;
Jmax = 200;
pop = Boltzmann(Temperature);               
delta_omega = delta_alpha * E0^2 / B/2;   
EstimateJmax = double(int16(pulse_FWHM / (1/(B0*1e2)/c / delta_omega / 2/pi))); % How to get the Jmax

[qq] = max(pop);
thre = qq*0.01;
ind = 1 : length(pop);
[qq2, Jhighpop] = max(ind(pop>thre));
Jmax = ceil((Jhighpop+EstimateJmax)/10) * 10;  % Set Jmax 10~20 higher than Estimated Jmax
upperJ = Jmax - 4;               % The highest rem. populated state due to temp
occ = zeros(upperJ, 2*upperJ+1, Jmax);  phase = occ;

%
global sigma; 
sigma1 = sqrt(pulse_FWHM^2 / (4*log(2)));  sigma = sigma1 * B / hbar;       
for Jstart = 0 : 1 : upperJ    
    Mdummy = 0;
    if mod(upperJ-Jstart, 10) == 0, disp(['Remaining J # is  ', num2str(upperJ-Jstart+1), '.']); end
    for Mstart = -Jstart : 1 : Jstart
        Mdummy = Mdummy + 1;
        c = coefficients(delta_omega, Jstart, Mstart);
        occ(Jstart+1, Mdummy, :) = abs(c);                
        phase(Jstart+1, Mdummy, :) = angle(c);       
    end
end

% Alignment Evolution over time of N2
t = 0.001 : 0.001 : 6.5;     % In unit of hbar/B
Jdummy = 0;  excitedJhighocc = 0; ciniJ = 0; ciniM = 0; 
for indexJtemp = 1 : upperJ   
    Jtemp = indexJtemp - 1;
    population = pop(indexJtemp);
    if population >= thre
        Mdummy = 0;
        for indexMtemp = 1 : (2*(Jtemp) + 1)   %this is an index for M
            M = indexMtemp - 1 - Jtemp;
            dummy = 0;            % Part for abs calculation
            for index = 1 : (Jmax-4)
                J = index - 1;
                if abs(occ(indexJtemp, indexMtemp, index)) > 0.05 && J > excitedJhighocc
                    excitedJhighocc = J;   ciniJ = Jtemp;   ciniM = M;
                end,  dummy = abs(occ(indexJtemp, indexMtemp, index)) .^ 2*cos2JM(J, M) + dummy;  end
            timeindep = dummy;             % Diagonal terms in sum over J&J'
            dummy = 0;            % Part for phase calculation
            for index = 1 : (Jmax-4)
                J = index - 1;
                phasep = phase(indexJtemp, indexMtemp, index) - phase(indexJtemp, indexMtemp, index+2);
                dummy = abs(occ(indexJtemp,indexMtemp,index)) * abs(occ(indexJtemp,indexMtemp,index+2))...
                    * (cos2Jp2M(J, M)) * cos((4*J+6).*t + phasep) + dummy;
            end,  timedep = dummy;       % Off diagonal terms in sum over J&J', only real part is considered
            Mdummy = 1/(2*Jtemp+1) * (2*timedep+timeindep) + Mdummy;  % Devide by M state number
        end
        Jdummy = population * Mdummy + Jdummy;
    end
end
cossquare = Jdummy;
[m1, p1] = max(cossquare);
[m2, p2] = min(cossquare);
% time = t * hbar / B * 1e12;  cossquare = Jdummy;
figure;  plot(t, cossquare, '-x');  grid on;  xlabel('time(ps)');  ylabel('<cos^2\theta>');
title([Mole, ' with E-density = ', num2str(E_density), 'J/cm^2 at Temperature = ', num2str(TemperatureK), 'K']);

%%  通过调整温度和激光强度，一个一个点算出来的，记录
Try_temp = [10, 20, 30, 40,43, 45, 50, 60, 70, 80, 90, 100, 110, 130];
try_tem_resul = [0.7340, 0.7449, 0.6942, 0.6498,0.6817 ,0.6768, 0.6645,...   % for tempretare=10:10:100K
                     0.6407, 0.6182, 0.5970, 0.5770, 0.5947, 0.5809, 0.5548];
Try_energy = [0.19612, 0.39224, 0.58836, 0.78448, 0.9806, 1.1767, 1.3728, 1.569, 1.7651, 1.9612, 2.1573,...
                      2.5496, 2.9418, 3.9224];
try_energy_resul  = [0.3406, 0.4917, 0.5481, 0.5724, 0.6638, 0.6798, 0.7304, 0.7449, 0.7427, 0.7835, 0.7838, ...
                        0.8095, 0.826, 0.8063];            % at tempretare = 20 k
figure; subplot(1,2,1);
plot(Try_temp, try_tem_resul);
subplot(1,2,2);            
plot(Try_energy, try_energy_resul);

%% 算角分布  角分布随时间的演化 图4.3(b)
resolution = [1, 41]; tic
Evolution = zeros(length(t), resolution(2));
for position2 = 1 : length(t)
    t2 = t(position2);
    wavepacketsq = 0;  wavepacket = 0;
    for indexJtemp = 1 : upperJ                                  %this is the Jtemp+1
        for indexMtemp = 1 : (2*(indexJtemp-1) + 1) %this is an index for 
            Jtemp = indexJtemp - 1;
            M = indexMtemp - 1 - Jtemp; 
            population = pop(indexJtemp) / (2*Jtemp+1);       

            wavepacket = 0;    
            for index = 1 : (Jmax)
                J = index-1;
                if occ(indexJtemp, indexMtemp, index) ~= 0
                    wavepacket = occ(indexJtemp,indexMtemp,index) * spharm(J,M,resolution,1)...
                                          * exp(-1i*J*(J+1)*t2 + 1i*phase(indexJtemp,indexMtemp,index)) + wavepacket;
                end
            end
            wavepacketsq = population * abs(wavepacket).^2 + wavepacketsq;
        end
    end
    distrib_evo = sum(wavepacketsq, 2);
    Evolution(position2, :) = distrib_evo; toc
end
% THETA = linspace(0, pi, resolution(2));
% [m1, p1] = max(cossquare);  [m2, p2] = min(cossquare);

%%    仅仅计算最大最小准直度时的角分布
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
THETA = linspace(0, pi, resolution(2)); distrib = sum(wavepacketsq, 2);
%  存文件单独再写一下
% filenameI = num2str(E_density);
% filenameT = num2str(TemperatureK);
% filename_cos = strcat('Cossquare-', Mole,'-Engery_density_', filenameI, '-Temp_', filenameT, 'K.dat');
% savevar(:, 1) = time;
% savevar(:, 2) = cossquare;
% save(filename_cos, 'savevar', '-ASCII');
% 
% filename_evol = strcat('Evolution-', Mole, '-Engery_density_', filenameI, '-Temp_', filenameT, 'K.dat');
% save(filename_evol, 'Evolution', '-ASCII');

%% figure 2.2 Alignment evolution of N2 over time
load('N2_alignment2.2.mat');
figure; set(gcf, 'position', [900, 90, 900, 700]);    
ax(11) = axes('Position', [0.08, 0.57, 0.37, 0.34]);
ax(12) = axes('Position', [0.57, 0.57, 0.38, 0.34]);
ax(2) = axes('Position', [0.1, 0.12, 0.85, 0.33]);

plot(ax(11), time(1:5800), cossquare(1:5800), '-');  grid(ax(11), 'on');  
xlabel(ax(11), 'time(ps)');  ylabel(ax(11), '<cos^2\theta>');  text(ax(11), 0.8, 0.47, '(A)', 'color', 'k');
title(ax(11), [Mole, ' with ', num2str(E_density), 'J/cm^2 at Temperature = ', num2str(TemperatureK), 'K']);
set(ax(11), 'FontSize', 14); 

% plot(ax(12), THETA, distrib, '-x'); grid(ax(12), 'on');   % 只画两张图时的角分布图
plot(ax(12), THETA, Evolution(p2, :), '-x'); grid(ax(12), 'on');
xlabel(ax(12), 'theta(θ)'); ylabel(ax(12), '归一化分布');  text(ax(12), 0.08, 0.11, '(B)', 'color', 'k');
title(ax(12), ['Max-cos^2  = ', num2str(m1), ', T = ', num2str(TemperatureK), 'K', ]);
set(ax(12), 'FontSize', 14); 

imagesc(ax(2), time, THETA, Evolution');    text(ax(2), 0.4, 0.2, '(C)', 'color', 'k');
xlabel(ax(2), 'Time delay (ps)'); ylabel(ax(2), 'The entire angle'); colorbar(ax(2));
set(gca, 'FontSize', 14);  

%% figure 2.2 Alignment evolution of  I2 over time
load('I2_alignment2.2.mat');
figure; set(gcf, 'position', [900, 90, 900, 400]);    
ax(11) = axes('Position', [0.08, 0.15, 0.37, 0.74]);
ax(12) = axes('Position', [0.57, 0.15, 0.38, 0.74]);
% ax(2) = axes('Position', [0.1, 0.12, 0.85, 0.33]);

plot(ax(11), time(1:5000), cossquare(1:5000), '-');  grid(ax(11), 'on');  
xlabel(ax(11), 'time(ps)');  ylabel(ax(11), '<cos^2\theta>');  text(ax(11), 50, 0.55, '(A)', 'color', 'k');
title(ax(11), [Mole1, ' with ', num2str(E_density), 'J/cm^2 at Temperature = ', num2str(TemperatureK), 'K']);
set(ax(11), 'FontSize', 14); 

plot(ax(12), THETA, distrib, '-x'); grid(ax(12), 'on');
% plot(ax(12), THETA, Evolution(p2, :), '-x'); grid(ax(12), 'on');
xlabel(ax(12), 'theta(θ)'); ylabel(ax(12), '归一化分布');  text(ax(12), 0.1, 0.14, '(B)', 'color', 'k');
title(ax(12), ['Max-cos^2  = ', num2str(m1), ', T = ', num2str(TemperatureK), 'K', ]);
set(ax(12), 'FontSize', 14); 

% imagesc(ax(2), time, THETA, Evolution');    text(ax(2), 0.4, 0.2, '(C)', 'color', 'k');
% xlabel(ax(2), 'Time delay (ps)'); ylabel(ax(2), 'The entire angle'); colorbar(ax(2));
% set(gca, 'FontSize', 14);  

%% figure 2.2 Alignment evolution of  CS2 over time
load('CS2_alignment2.2.mat');
figure; set(gcf, 'position', [900, 90, 900, 400]);    
ax(11) = axes('Position', [0.08, 0.15, 0.37, 0.74]);
ax(12) = axes('Position', [0.57, 0.15, 0.38, 0.74]);
% ax(2) = axes('Position', [0.1, 0.12, 0.85, 0.33]);

plot(ax(11), time(1:5000), cossquare(1:5000), '-');  grid(ax(11), 'on');  
xlabel(ax(11), 'time(ps)');  ylabel(ax(11), '<cos^2\theta>');  text(ax(11), 50, 0.55, '(A)', 'color', 'k');
title(ax(11), [Mole2, ' with ', num2str(E_density), 'J/cm^2 at Temperature = ', num2str(TemperatureK), 'K']);
set(ax(11), 'FontSize', 14); 

plot(ax(12), THETA, distrib, '-x'); grid(ax(12), 'on');
% plot(ax(12), THETA, Evolution(p2, :), '-x'); grid(ax(12), 'on');
xlabel(ax(12), 'theta(θ)'); ylabel(ax(12), '归一化分布');  text(ax(12), 0.1, 0.14, '(B)', 'color', 'k');
title(ax(12), ['Max-cos^2  = ', num2str(m1), ', T = ', num2str(TemperatureK), 'K', ]);
set(ax(12), 'FontSize', 14); 

% imagesc(ax(2), time, THETA, Evolution');    text(ax(2), 0.4, 0.2, '(C)', 'color', 'k');
% xlabel(ax(2), 'Time delay (ps)'); ylabel(ax(2), 'The entire angle'); colorbar(ax(2));
% set(gca, 'FontSize', 14);  
%%  准直影响因素画图
Try_temp = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110];
try_tem_resul = [0.7340, 0.7449, 0.6942, 0.6498,0.6645,...   % for temperature=10:10:100K
                     0.6407, 0.6182, 0.5970, 0.5770, 0.5947, 0.5809];
Try_energy = [0.19612, 0.39224, 0.58836, 0.78448, 0.9806, 1.1767, 1.3728, 1.569, 1.7651, 1.9612, 2.1573,...
                      2.5496, 2.9418, 3.9224];
try_energy_resul  = [0.3406, 0.4917, 0.5481, 0.5724, 0.6638, 0.6798, 0.7304, 0.7449, 0.7427, 0.7835, 0.7838, ...
                        0.8095, 0.826, 0.8063];            % at temperature = 20 k
figure;  subplot(1,2,1);
plot(Try_temp, try_tem_resul, '-x'); grid on
xlabel('T(K)'); ylabel('Max <cos^2\theta>'); text(10, 0.735, '(A)', 'color', 'k');
title('Max- <cos^2\theta> with temperature ');set(gca, 'FontSize', 13); 

subplot(1,2,2);            
plot(Try_energy, try_energy_resul, '-x'); grid on 
xlabel('Energy(J/cm^2)'); ylabel('Max <cos^2\theta>');
title('Max- <cos^2\theta> with Energy'); text(0.3, 0.86, '(B)', 'color', 'k');set(gca, 'FontSize', 13); 
