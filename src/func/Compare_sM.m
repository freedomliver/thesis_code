function [Exp_sM, Factor] = Compare_sM(s,EffExpI,SimIA,SimIM,RadiusRange)
% comparing experimental sM and simulative sM
% 2021/03/10 qifengfeng
%
% s: the s range; 1D array
% EffExpI: experimental intensity; 1D array
% SimIA: simulation of IA; 1D array
% SimIM: simultation of IM; 1D array
% RadiusRange: the range of effective pixel containing in one diffraction
% pattern; 1D array
%
% Exp_sM: experimental sM; 1D array

%% finding where simulative IM is equal to zero
Multiplier=s.*0;
Multiplier(2:end)=SimIM(1:end-1);
FindZeroIM=SimIM.*Multiplier; % Multiplication of adjacent data
ZeroPos=find(FindZeroIM<0); % where IM is equal to zero

%% fitting IA in experiment (using (IA+IM)/IA) best way

ExpData=EffExpI./SimIA; % fitting data; /SimIA is to avoid rapid attenuation and for better fitting
% fitting points
FitData=ExpData(ZeroPos);
Fit_s=s(ZeroPos);
% polyfit
Fit_coff=polyfit(Fit_s,FitData,length(ZeroPos)-1);
% fitting result
FitResult=polyval(Fit_coff,s);  % Fit IA

%% calculating sM in experiment

Exp_sM=(ExpData-FitResult).*s; % sM in experiment
Sim_sM=s.*SimIM./SimIA; % sM of simulation
% normalization
Factor=Sim_sM(RadiusRange(1))/Exp_sM(RadiusRange(1))*0.5;
Exp_sM=Factor.*Exp_sM;

figure
plot(s(RadiusRange),Exp_sM(RadiusRange),'b-','LineWidth',1.8)
hold on
plot(s,Sim_sM,'r-','LineWidth',1.8)
xlabel('s ($\AA^{-1}$)','Interpreter','latex')
ylabel('sM (arb. unit)')
xlim([0,10])
grid on
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',30)
set(gcf,'unit','normalized','position',[0.3,0.3,0.45,0.55]);
set (gca,'position',[0.2,0.2,0.64,0.64] );
head=legend('experiment','simulation');
set(head,'FontName','Times New Roman','FontSize',20,'FontWeight','normal')

end