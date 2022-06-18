function [AveragePattern,AverageMask,StdPattern]=AveragePattern(Pattern,Mask,NumCycles,ThresholdStd)
% calculating the mean pattern of the same delay and removing the hot pixels
% 2020/10/20 qifengfeng
% 
% Pattern: patterns with same time delay; N*N*number array
% Mask: mask patterns with same time delay; N*N*number array
% NumCycles: the number of average cycles to remove the signal spikes from random events
% ThresholdStd: signals within ThresholdStd standard deviations are reserved
%
% AveragePattern: the average pattern of all patterns with same time delay
% AverageMask: the mask pattern of the average pattern
% StdPattern: the standard deviations of each pixel

for ii=1:NumCycles
    
    % calculating the average value of each pixel
    AveragePattern=sum(Pattern.*Mask,3)./sum(Mask,3);  % result contianing NaN
    % calculating the standard deviation of each pixel
    StdPattern=sqrt(sum(((Pattern-AveragePattern).*Mask).^2,3)./sum(Mask,3));
    % calculating the upper limit and lower limit of valid data
    UpperLimit=AveragePattern+StdPattern.*ThresholdStd;
    LowerLimit=AveragePattern-StdPattern.*ThresholdStd;
    % find valid data
    ValidData=(Pattern<=UpperLimit & Pattern>=LowerLimit);
    % refresh the mask pattern
    Mask = Mask.*ValidData;
    
end

% calculating the average value of each pixel
AveragePattern=sum(Pattern.*Mask,3)./sum(Mask,3);
% calculating the standard deviation of each pixel
StdPattern=sqrt(sum(((Pattern-AveragePattern).*Mask).^2,3)./sum(Mask,3));

% valid datas of each pixel
AverageMask=sum(Mask,3);
% for each pixel, the number of valid datas take more than a third to be effective
SizePattern=size(Pattern);
if length(SizePattern)==2
    SizePattern(3)=1;
end
AverageMask=double(AverageMask>(SizePattern(3)/30));

% updata average pattern and std pattern
AveragePattern=AveragePattern./(AverageMask./AverageMask); % avoid Inf
StdPattern=StdPattern./(AverageMask./AverageMask);

end