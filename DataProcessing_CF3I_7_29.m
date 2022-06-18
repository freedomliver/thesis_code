% Data processing (Gas UED) for 2021-7/29-30
clc;  clear;  addpath(genpath('D:\Program Files\MATLAB\cal\Diffraction_qi\src'));  % 加载程序包位置
MoleculeFile='D:\Program Files\MATLAB\cal\Diffraction_qi\src\molecule\CF3I.xyz';
RadiusRange = 1 : 450;  sPixel = 0.0264;        % s value of each pixel (angstrom^-1)
s = RadiusRange.*sPixel - sPixel/2; 
[SimIA, SimIM] = SimIAIM(MoleculeFile, s);  % simultion of IA and IM

% RadiusRange = 1 : 4500;  sPixel = 0.0264;        % s value of each pixel (angstrom^-1)
% s = RadiusRange.*sPixel / 10 - sPixel/2; 

% BackPath = uigetdir('F:\');                                % loading background
% BackPathListing = dir([BackPath, '\*.tiff']);
% BackPatterns = zeros(1024, 1024, numFolder);
% BackpatternNames(:) = strcat({BackPathListing.folder}, '\', {BackPathListing.name}); 
% for ii = 1 : length(BackpatternNames)  
%     BackPatterns(:, :, ii) = imread(BackpatternNames{ii});  end


% BackFile = 'D:\CF3I_Data\2021-07-29 CF3I\back.spe';       % loading background  提前设定大小试一下。
% BackPatterns = ReadSPE(BackFile);


% save('loading.mat');

%% loading data  加载数据包
clc; clear; 
% load('loading.mat');

% start_num = 122;  finish_num = 200;                 % Begin and the end of the sub folder
% [path, filefolders, numFolder, numPattern] = LoadingFolder('D:\CF3I_Data\', start_num);%, finish_num);
[path, filefolders, numFolder, numPattern] = LoadingFolder('D:\CF3I_Data\');
[PatternNames, DelayTime] = LoadingName(path, filefolders, numFolder, numPattern);
TimeDelay = mean(DelayTime);                                           % calculating time delay
TimeDelay=(TimeDelay-TimeDelay(1)) ./ 1e7 ./ 3e8 .* 1e15;

%% find hole center  寻找衍射图的中心位置，中心圆圈只和荧光屏有关系，不会随着电子抖动发生改变。
SinglePattern = double(imread(PatternNames{1,1}));
[numRowPix, numColPix] = size(SinglePattern);   
SizePix = [numRowPix, numColPix];                                     % 用一个数据获得衍射图大小

HolePattern = zeros(numRowPix, numColPix, numPattern);% 每个文件夹第一张图
for ii = 1 : numPattern,  HolePattern(:, :, ii) = double(imread(PatternNames{ii, 1}));  end 
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

%%  画单幅与160幅的对比图
figure;
ax(1) = axes('Position',[0.08, 0.1, 0.4, 0.8]);
pcolor(SinglePattern),  shading interp,  colormap(jet), colorbar;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
title('Just one pattern signal');  text(100, 900, '(A)', 'Color','w');

ax(2) = axes('Position',[0.55, 0.1, 0.4, 0.8]);
pcolor(median(HolePattern, 3)),  shading interp,  colormap(jet),colorbar;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
title('Accumulated pattern');text(100, 900, '(B)', 'Color','w');

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
numPix = 851;  Range = (numPix-1)/2; 
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
%         CenterCol = 443;
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

%% Averaging back pattern
numPix = 851;  Range = (numPix-1)/2; 
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

%% Static pattern before time zero
BeforeZero = 1:6 ; % static patterns index
[StaticPattern, StaticMask, ~] = AveragePattern(TimePattern(:,:,BeforeZero), TimeMask(:,:,BeforeZero), 1, 3);

TimeDiff = TimePattern;  % Calculating the difference patterns (ΔIM)
for ii = 1 : numPattern
    tmpDiff = TimePattern(:,:,ii) - StaticPattern;
    tmpDiff((TimeMask(:,:,ii) .* StaticMask) == 0) = NaN;
    TimeDiff(:,:,ii) = tmpDiff;
end

% figure
% pcolor(TimeDiff(:,:,end)), shading interp, colormap(jet), caxis([-25,25])
% 


%%
  addpath(genpath('D:\Program Files\MATLAB\cal\Diffraction_qi\src'));  % 加载程序包位置
tttt = TimeDiff;
% TimeDiff_para = TimeDiff(380: 478, :,:);
TimeDiff_para = TimeDiff(395: 463, :,:);
TimeDiff_vertical = TimeDiff(:,395: 463,:);
% figure;               
% pcolor(TimeDiff_para(:, :, 1)),  shading interp,  colormap(jet), caxis([-25,25])
%% Radial average
 TimeDiff = TimeDiff_vertical;
[~, numIM, ~, ~] = RadialAverage(TimeDiff(:,:,1),3);
numRadius = size(numIM, 2);

DeltaIM = zeros(numPattern, numRadius);  numIM = DeltaIM;  StdIM = DeltaIM;  % initialization 
ModifiedDiff = TimeDiff;

CorrectFactor = 3; % Any pixel that are more than 3 stds away from the mean is removed
for ii = 1 : numPattern
    [DeltaIM(ii,:), numIM(ii,:), StdIM(ii,:), ModifiedDiff(:,:,ii)] = RadialAverage(TimeDiff(:,:,ii), CorrectFactor);
end
% 
% figure  % test ModifiedDiff
% surf(ModifiedDiff(:,:,end)), shading interp, colormap(jet)
% disp('===Difference patterns with different time delay have been acquired===')

%% acquiring effective s range
% RadiusRange = 1 : 420;  sPixel = 0.0264;        % s value of each pixel (angstrom^-1)
% s = RadiusRange.*sPixel - sPixel/2; 
% [SimIA, SimIM] = SimIAIM(MoleculeFile, s);  % simultion of IA and IM
EffDeltaIM = DeltaIM(:, RadiusRange);    % effective pixels
EffnumIM = numIM(:, RadiusRange);
EffStdIM = StdIM(:, RadiusRange);

[StaticI, ~, ~, ~] = RadialAverage(StaticPattern, CorrectFactor);  % Percent Difference
figure
plot(s, StaticI(1, RadiusRange))
PD = DeltaIM(:, RadiusRange) ./ StaticI(1, RadiusRange) .* 100;

figure
plot(s, PD(1, RadiusRange), 'LineWidth', 1.5)
hold on
for ii = 10:13 ,  plot(s,PD(ii,RadiusRange),'LineWidth',1.5);  end; hold on
%  plot(s, PD(9,RadiusRange),'LineWidth',1.5); 
xlabel('s ($\AA^{-1}$)', 'Interpreter', 'latex');  ylabel('PD (%)');  xlim([1,10]);
grid on, box on
set(gca, 'linewidth', 1.2);  set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
set(gcf, 'unit', 'normalized', 'position', [0.3, 0.3, 0.45, 0.55]);
set(gca, 'position', [0.2, 0.2, 0.64, 0.64]);  legend

% figure
% plot(s, PD(1, RadiusRange), 'LineWidth', 1.5);
% xlabel('s ($\AA^{-1}$)', 'Interpreter', 'latex');  ylabel('PD (%)');  xlim([1,8]);
% grid on, box on
% set(gca, 'linewidth', 1.2);  set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
% set(gcf, 'unit', 'normalized', 'position', [0.3, 0.3, 0.45, 0.55]);
% set(gca, 'position', [0.2, 0.2, 0.64, 0.64]);  legend

%% 查看信号变化强度
PD1 = mean(PD(1 : 2, RadiusRange));     % 
PD2 = mean(PD(4:6, RadiusRange));
figure
plot(s, PD1, 'LineWidth', 1.5);  %
hold on;  plot(s, PD2, 'LineWidth', 1.5);
xlabel('s ($\AA^{-1}$)','Interpreter','latex');  ylabel('PD (%)');  xlim([1,6]);
grid on, box on
set(gca, 'linewidth', 1.2);  set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
set(gcf, 'unit', 'normalized', 'position', [0.3, 0.3, 0.45, 0.55]);
set(gca, 'position', [0.2, 0.2, 0.64, 0.64]);  legend

%% Ensure the diffraciton pattern before time zero
sMRadiusRange = 60 : 410;
Exp_sM = Compare_sM(s, StaticI(RadiusRange), SimIA, SimIM, sMRadiusRange);

% calculating experimental ΔsM and ΔPDF
ExpDeltaSM = EffDeltaIM ./ SimIA .* s;   % ΔsM
r = 0 : 0.001 : 10; % angstrom 
ExpDeltaPDF = 0 .* r;       % ΔPDF
ExpDeltaPDF = repmat(ExpDeltaPDF, numPattern, 1);
k = 0.04;
for ii = 1 : length(r)
    tmp_r = r(ii);    % the data of NaN is like the data 0
    ExpDeltaPDF(:, ii) = nansum(ExpDeltaSM .* sin(s.*tmp_r) .* exp(-1.*k.*(s.^2)) .* (s(2)-s(1)), 2); 
end
disp('===Experimental ΔsM and ΔPDF have been acquired===')

% figure of experimental ΔsM and ΔPDF
% figure
% pcolor(s, TimeDelay.*2 - 160, ExpDeltaSM),  shading flat,  colormap(jet)
% xlabel('s ($\AA^{-1}$)', 'Interpreter', 'latex');  ylabel('Delay (fs)');  xlim([0, 8]);
% grid on, box on
% set(gca, 'linewidth', 1.2);  set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
% set(gcf, 'unit', 'normalized', 'position', [0.3, 0.3, 0.45, 0.55]);
% set(gca, 'position', [0.2, 0.2, 0.64, 0.64]);  legend
ExpDeltaPDF_para = ExpDeltaPDF;
% ExpDeltaPDF_perp = ExpDeltaPDF;


fig = figure; set(gcf, 'position', [600, 90, 1000, 400]);    % subplot[5,2]. IT, IA, IM, sM and f(r) for each side line
ax(1) = axes('Position', [0.1, 0.12, 0.4, 0.8]);
ax(2) = axes('Position', [0.55, 0.12, 0.4, 0.8]);

pcolor(ax(1), r, TimeDelay.*2 - 160, ExpDeltaPDF_para),  shading(ax(1),'flat'),  colormap(ax(1),'jet')
xlabel(ax(1), 'r($\AA$)', 'Interpreter', 'latex');  ylabel(ax(1), 'Delay (fs)');  colorbar(ax(1)); 
xlim(ax(1),[1,6]); title(ax(1),'The parallel pdf signal');

pcolor(ax(2), r, TimeDelay.*2 - 160, ExpDeltaPDF_perp),  shading flat,  colormap(jet)
xlabel(ax(2),'r($\AA$)', 'Interpreter', 'latex');  ylabel(ax(2),'Delay (fs)'); colorbar(ax(2))
xlim([1,6]); title(ax(2),'The perpendicular pdf signal');

% 通过将其他区域都删掉，只留竖直和水平的信号看看有没有各向异性。

%% results of experimental ΔsM and ΔPDF
figure   % ΔsM
plot(s, ExpDeltaSM(1,:), 'LineWidth', 1.5);  hold on
for ii = 2:2 : 6,  plot(s, ExpDeltaSM(ii,:), 'LineWidth', 1.5);  end
% for ii = 30 : 35,  plot(s, ExpDeltaSM(ii,:), '--', 'LineWidth', 1.5);  end
xlabel('s ($\AA^{-1}$)', 'Interpreter', 'latex');  ylabel('\DeltasM (arb. unit)');  xlim([1,10]);
grid on, box on;
set(gca, 'linewidth', 1.2);  set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
set(gcf, 'unit', 'normalized', 'position', [0.3, 0.3, 0.45, 0.55]);
set(gca, 'position', [0.2, 0.2, 0.64, 0.64]);  legend


figure  % ΔPDF
plot(r, ExpDeltaPDF(2,:), 'LineWidth', 1.5);  hold on
for ii = 19:22,  plot(r, ExpDeltaPDF(ii,:), 'LineWidth', 1.5);  end
% for ii = 18:2 : 24,  plot(r, ExpDeltaPDF(ii,:), '--', 'LineWidth', 1.5);  end
xlabel('r ($\AA$)', 'Interpreter', 'latex');  ylabel('ΔPDF (arb. unit)');  xlim([1, 6]);
grid on, box on
set(gca, 'linewidth', 1.2);  set(gca, 'FontName', 'Times New Roman', 'FontSize', 30);
set(gcf, 'unit', 'normalized', 'position', [0.3, 0.3, 0.45, 0.55]);
set(gca, 'position', [0.2, 0.2, 0.64, 0.64]);  legend

% save('7_30_Scan12_wu84-106_vertical.mat');

%% 

figure;
pcolor(SinglePattern),  shading interp,  colormap(jet),  colorbar
xlabel('Pixels'), ylabel('Pixels');  title('Just one Pattern signal');
set(gca,'FontSize', 19);
set(gcf, 'position', [100, 100, 800, 800]);  

pcolor(ax(1), SinglePattern);  colorbar(ax(1));  %text(ax(1), -9, -9, '(A)', 'color', 'r');
xlabel(ax(1), 'sx ($\AA^{-1}$)','Interpreter','latex'), ylabel(ax(1), 'sy ($\AA^{-1}$)','Interpreter','latex')
set(gca, 'xtick', -10 : 5 : 10);  set(gca, 'xticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'ytick', -10 : 5 : 10);  set(gca, 'yticklabel', [10, 5, 0, 5, 10]); 
set(gca, 'FontSize', 8);  title(ax(1), 'Perfect Alignment-IT Pattern');

figure;
pcolor(median(HolePattern, 3)),  shading interp,  colormap(jet),  colorbar
xlabel('Pixels'), ylabel('Pixels');  title('Accumulate pattern');
set( gca,'FontSize', 19);
set(gcf, 'position', [100, 100, 800, 800]);