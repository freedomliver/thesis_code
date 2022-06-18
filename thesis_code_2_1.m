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

%% figure 2.1 perfect align and random orientation
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

%%
save('article_chapter2.1.mat');
load('article_chapter2.1.mat');

%%  处理和模拟碘分子静态数据
load('I2_static2.1.mat');
%%
% clc;  clear; 
addpath(genpath('D:\Program Files\MATLAB\cal\src')); % add path for program package

% loading file
TotalFile='D:\气态\I2-10s-45c.spe';  % molecule+background signal
BackFile='D:\气态\I2-10s-45c-back.spe'; % background signal 
TotalIMG=ReadSPE(TotalFile);
BackIMG=ReadSPE(BackFile);

disp('===Pattern files have been loaded===')

%% finding hole center, pattern center and creating mask pattern
% finding hole
HoleFinderPattern = median(TotalIMG,3);
HoleStrengthMax = 700; % the upper limit of the pattern strength in the hole; for converting the pattern to a binary image
HoleDeviation = 60; % HoleDeviation: the rough pixels of the hole edge off the center of the pattern
[HoleCenterRow,HoleCenterCol,HoleRadius] = FindHoleCenter(HoleFinderPattern,HoleStrengthMax,HoleDeviation);

% single mask pattern; mask data in the hole is 0
[numRowPix,numColPix,~] = size(TotalIMG);
MaskPatternSingle=ones(numRowPix,numColPix);
for ii=1:numRowPix
    for jj=1:numColPix
        if sqrt((ii-HoleCenterRow)^2+(jj-HoleCenterCol)^2)<=(HoleRadius+1)
            MaskPatternSingle(ii,jj)=0;
        end
    end
end
% MaskPatternSingle(866,839)=0; % this pixel in EMCCD has been broken (row866,column839)

% remove the region contianing dark current
DarkCurrentFig = figure;
pcolor(HoleFinderPattern), shading interp, colormap(jet), caxis([450,700]) 
set(gca,'FontName','Times New Roman','FontSize',12)
set(gcf,'position',[100,100,800,800]);
axis equal

Flag=0;
Choice=questdlg('Do you want to remove dark current', 'Removing Dark Current','Yes','No','Yes');
if strcmp(Choice,'Yes') == 1
    Flag=1;
end
    
while Flag==1
    DarkCurrent=getrect(DarkCurrentFig);
    DarkColumn=round(DarkCurrent(1)):round(DarkCurrent(1))+round(DarkCurrent(3));
    DarkRow=round(DarkCurrent(2)):round(DarkCurrent(2))+round(DarkCurrent(4));
    MaskPatternSingle(DarkRow,DarkColumn)=0;
    
    Choice=questdlg('Do you want to continue', 'Removing Dark Current','Yes','No','No');
    if strcmp(Choice,'No') == 1
        Flag=0;
    end     
end

% finding the center of diffraction pattern
[CenterCol,CenterRow]=getpts(DarkCurrentFig);
CenterCol=round(CenterCol);
CenterRow=round(CenterRow);

disp('===pattern center has been found and mask pattern has been created===')

%% cutting patterns
% determining the region
numPix=901; % modifying patterns to NumPix*NumPix; odd number
Range=(numPix-1)/2; % center+-Range
RowRange=(CenterRow-Range):(CenterRow+Range);
ColRange=(CenterCol-Range):(CenterCol+Range);

% modifying patterns
% the first pattern in .spe file collected by EMCCD is not accurate
TotalPattern=TotalIMG(RowRange,ColRange,2:end);
BackPattern=BackIMG(RowRange,ColRange,2:end);
MaskSingle=MaskPatternSingle(RowRange,ColRange);

disp('===Patterns have been aligned center===')

%% averging patterns

TotalSzie=size(TotalPattern);
Mask=repmat(MaskSingle,1,1,TotalSzie(3));
% calculating the mean pattern of the same delay and removing the hot pixels
numCycles=3; % the number of average cycles to remove the signal spikes from random events
ThresholdStd=3; % signals within ThresholdStd standard deviations are reserved
% TimePattern and TimeStd contain NaN
[TotalAverage,TotalMask,TotalStd]=AveragePattern(TotalPattern,Mask,numCycles,ThresholdStd);

% test the quality of the average
figure()
subplot(1,3,1)
pcolor(TotalAverage), shading interp, colormap(jet)
subplot(1,3,2)
pcolor(TotalMask), shading interp, colormap(jet)
subplot(1,3,3)
surf(TotalStd), shading interp, colormap(jet)
set(gcf,'unit','centimeters','position',[1,10,50,12])

BackSzie=size(BackPattern);
Mask=repmat(MaskSingle,1,1,BackSzie(3));
% calculating the mean pattern of the same delay and removing the hot pixels
% TimePattern and TimeStd contain NaN
[BackAverage,BackMask,BackStd]=AveragePattern(BackPattern,Mask,numCycles,ThresholdStd);

disp('===Average patterns have been acquired===')

%% acquiring IA+IM in experiment

ExpPattern = TotalAverage-BackAverage; % IA+IM; molecular pattern
% pcolor(ExpPattern), shading interp, colormap(jet)
% set(gca,'FontName','Times New Roman','FontSize',12)
% set(gcf,'position',[100,100,800,800]);

% median filter 
ExpPattern= medfilt2(ExpPattern,[9,9]);
pcolor(ExpPattern), shading interp, colormap(jet)

% radial average
correct=3;
[ExpI,ExpNumI,ExpPatternStd,ExpPatternModified]=RadialAverage(ExpPattern,correct); % IA+IM (existing NAN)

RadiusRange=1:460; % the length of effective pixel containing in one diffraction pattern
sPixel=0.025; % s value of each pixel (angstrom^-1)
s=RadiusRange.*sPixel-sPixel/2; % s

% effective pixels
EffExpI=ExpI(RadiusRange);
figure;
plot(s,EffExpI)

disp('===Radial average has been acquired===')

%% acquiring scattering amplitude
MoleculeFile='D:\Program Files\MATLAB\cal\src\molecule\I2.xyz';
[SimIA,SimIM]=SimIAIM(MoleculeFile,s);

% figure;
% plot(s,SimIA); hold on;
% plot(s,SimIM); legend
%% finding where simulative IM is equal to zero

Multiplier=s.*0;
Multiplier(2:end)=SimIM(1:end-1);
FindZeroIM=SimIM.*Multiplier; % Multiplication of adjacent data
ZeroPos=find(FindZeroIM<0); % where IM is equal to zero

disp('===IM of simulation equal to zero have been found===')

%% fitting IA in experiment (using (IA+IM)/IA) best way

ExpData=EffExpI./SimIA./1e16; % fitting data; /SimIA is to avoid rapid attenuation and for better fitting
% fitting points
FitData=ExpData(ZeroPos);
Fit_s=s(ZeroPos);
% polyfit
Fit_coff=polyfit(Fit_s,FitData,length(ZeroPos)-1);
% fitting result
FitResult=polyval(Fit_coff,s); 

figure
plot(s,ExpData,'b-','Linewidth',1.8)
hold on
plot(Fit_s,FitData,'mo','MarkerSize',5,'LineWidth',1.5)
plot(s,FitResult,'r','Linewidth',1.8)
xlabel('s ($\AA^{-1}$)','Interpreter','latex')
ylabel('R (arb. unit)')
xlim([0,13])
grid on
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',30)
set(gcf,'unit','normalized','position',[0.3,0.3,0.45,0.55]);
set (gca,'position',[0.2,0.2,0.64,0.64] );
head=legend('I','Fit points','Fit background');
set(head,'FontName','Times New Roman','FontSize',20,'FontWeight','normal');legend

disp('===IA has been fitted===')

 %% fitting IA in experiment (using IA+IM)

ExpData=EffExpI./1e3; % fitting data; /SimIA is to avoid rapid attenuation and for better fitting
% fitting points
FitData=ExpData(ZeroPos);
Fit_s=s(ZeroPos);
% polyfit
Fit_coff=polyfit(Fit_s,FitData,length(ZeroPos)-1);
% fitting result
FitResult=polyval(Fit_coff,s); 

figure
plot(s,ExpData,'b-','Linewidth',1.8)
hold on
plot(Fit_s,FitData,'mo','MarkerSize',5,'LineWidth',1.5)
plot(s,FitResult,'r','Linewidth',1.8)
xlabel('s ($\AA^{-1}$)','Interpreter','latex')
ylabel('I (arb. unit)')
xlim([0,13])
grid on
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',30)
set(gcf,'unit','normalized','position',[0.3,0.3,0.45,0.55]);
set (gca,'position',[0.2,0.2,0.64,0.64] );
head=legend('I','Fit points','Fit background');
set(head,'FontName','Times New Roman','FontSize',20,'FontWeight','normal')

disp('===IA has been fitted===')

%% fitting IA in experiment (using (IA+IM)/IA*s)

ExpData=EffExpI./SimIA./1e16.*s; % fitting data; /SimIA is to avoid rapid attenuation and for better fitting
% fitting points
FitData=ExpData(ZeroPos);
Fit_s=s(ZeroPos);
% polyfit
Fit_coff=polyfit(Fit_s,FitData,length(ZeroPos)-1);
% fitting result
FitResult=polyval(Fit_coff,s); 

figure
plot(s,ExpData,'b-','Linewidth',1.8)
hold on
plot(Fit_s,FitData,'mo','MarkerSize',5,'LineWidth',1.5)
plot(s,FitResult,'r','Linewidth',1.8)
xlabel('s ($\AA^{-1}$)','Interpreter','latex')
ylabel('R (arb. unit)')
xlim([0,13])
grid on
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',30)
set(gcf,'unit','normalized','position',[0.3,0.3,0.45,0.55]);
set (gca,'position',[0.2,0.2,0.64,0.64] );
head=legend('I','Fit points','Fit background');
set(head,'FontName','Times New Roman','FontSize',20,'FontWeight','normal')

disp('===IA has been fitted===')

%% calculating sM in experiment
RadiusRange=50:400; % the range of effective pixel containing in one diffraction pattern

Exp_sM=(ExpData-FitResult).*s; % sM in experiment
Sim_sM=s.*SimIM./SimIA; % sM of simulation
% normalization
Exp_sM=Sim_sM(RadiusRange(1))/Exp_sM(RadiusRange(1)).*Exp_sM.*1;

figure
plot(s(RadiusRange),Exp_sM(RadiusRange),'b-','LineWidth',1.8)
hold on
plot(s,Sim_sM*9,'r-','LineWidth',1.8)
xlabel('s ($\AA^{-1}$)','Interpreter','latex')
ylabel('sM (arb. unit)')
xlim([0,13])
grid on
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',30)
set(gcf,'unit','normalized','position',[0.3,0.3,0.45,0.55]);
set (gca,'position',[0.2,0.2,0.64,0.64] );
head=legend('experiment','simulation');
set(head,'FontName','Times New Roman','FontSize',20,'FontWeight','normal')

%% calculating PDF
% Exp_sM of s<smin is replaced by Sim_sM
ModifiedFilt_sM=Exp_sM;
ModifiedFilt_sM(1:RadiusRange(1))=Sim_sM(1:RadiusRange(1));

figure
plot(s,ModifiedFilt_sM,'b-','LineWidth',1.8)
xlabel('s ($\AA^{-1}$)','Interpreter','latex')
ylabel('sM (arb. unit)')
xlim([0,13])
grid on
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',30)
set(gcf,'unit','normalized','position',[0.3,0.3,0.45,0.55]);
set (gca,'position',[0.2,0.2,0.64,0.64] );

% ΔPDF
r=0:0.001:10; % angstrom
ExpPDF=0.*r;
SimPDF=0.*r;
k=0.06;
for ii=1:length(r)
    ExpPDF(ii)=nansum(ModifiedFilt_sM.*sin(s.*r(ii)).*exp(-1.*k.*(s.^2)).*(s(2)-s(1))); 
    SimPDF(ii)=nansum(Sim_sM.*sin(s.*r(ii)).*exp(-1.*k.*(s.^2)).*(s(2)-s(1)));
end

figure
plot(r,ExpPDF,'b-','LineWidth',1.8); hold on
plot(r,SimPDF*4.5,'r-','LineWidth',1.8); grid on
xlabel('r ($\AA$)','Interpreter','latex'); ylabel('PDF (arb. unit)')
xlim([0,13])
grid on
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',30)
set(gcf,'unit','normalized','position',[0.3,0.3,0.45,0.55]);
set (gca,'position',[0.2,0.2,0.64,0.64] );
head=legend('experiment','simulation');
set(head,'FontName','Times New Roman','FontSize',20,'FontWeight','normal')

%%  画图的axis可以放在每一张图的前面，可以省掉标图轴
fig = figure; set(gcf, 'position', [700, 90, 800, 900]);    
ax(21) = axes('Position', [0.08, 0.57, 0.4, 0.35]);
ax(22) = axes('Position', [0.55, 0.57, 0.4, 0.35]);
sx_pat = -450:450;
imagesc(ax(21), sx_pat*0.0264, sx_pat*0.0264, TotalAverage), colorbar(ax(21)); 
xlabel(ax(21), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'), ylabel(ax(21), 'sy ($\AA^{-1}$)','Interpreter', 'latex');
title(ax(21), 'Experiment total pattern');
set(ax(21), 'FontSize', 13);

FitD_fig=SimIM(ZeroPos) * 7;
plot( s,EffExpI,'b-','Linewidth',1.8); hold on 
plot(s,SimIM*7); plot(s,SimIA*7); 
plot(Fit_s, FitD_fig,'mo','MarkerSize',5,'LineWidth',1.5);
legend('Experiment-data','Simulation-IM','Simulation-IA','Fitted zero-point');
xlabel(ax(22), 'sx ($\AA^{-1}$)', 'Interpreter', 'latex'),% ylabel(ax(22), 'sy ($\AA^{-1}$)', 'Interpreter', 'latex');
title(ax(22), 'Experiment and simulation result'); grid on;
set(ax(22), 'FontSize', 13);

ax(11) = axes('Position', [0.08, 0.08, 0.4, 0.38]);
plot(s(RadiusRange),Exp_sM(RadiusRange),'b-','LineWidth',1.8); hold on
plot(s,Sim_sM*9,'r-','LineWidth',1.8);xlabel('s ($\AA^{-1}$)','Interpreter','latex');
ylabel('sM (arb. unit)');grid on 
title('sM of experiment and simulation'); set(ax(11), 'FontSize', 13);
head=legend('experiment','simulation');

ax(12) = axes('Position', [0.55, 0.08, 0.4, 0.38]);
plot(r,ExpPDF,'b-','LineWidth',1.8); hold on
plot(r,SimPDF*4.5,'r-','LineWidth',1.8); grid on
xlabel('r ($\AA$)','Interpreter','latex'); ylabel('PDF (arb. unit)');
legend('experiment','simulation');
title('PDF of experiment and simulation');
set(ax(12), 'FontSize', 13);

text(ax(21), -10, -10, '(A)', 'color', 'k');
text(ax(22), 1.3, 700, '(B)', 'color', 'k');
text(ax(11), 0.6, 12.5, '(C)', 'color', 'k');
text(ax(12), 0.6, 3.5, '(D)', 'color', 'k');

