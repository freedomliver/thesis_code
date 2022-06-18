function [SimIA,SimIM] = SimIAIM(MoleculeFile,s)
% Simulation of IA and IM
% 2021/03/08 qifengfeng
%
% MoleculeFile: the .xyz file of the structure of molecule
% s: the s range; 1D array
%
% SimIA: simulation result of IA; 1D array
% SimIM: simulation result of IM; 1D array

%% scattering amplitude
DPWAfile='D:\Program Files\MATLAB\cal\src\DPWA\'; % the path of the file DPWA
[fH,~]=DPWAcpu(s,1,3000,DPWAfile);
[fHe,~]=DPWAcpu(s,2,3000,DPWAfile);
[fC,~]=DPWAcpu(s,6,3000,DPWAfile);
[fN,~]=DPWAcpu(s,7,3000,DPWAfile);
[fO,~]=DPWAcpu(s,8,3000,DPWAfile);
[fF,~]=DPWAcpu(s,9,3000,DPWAfile);
[fS,~]=DPWAcpu(s,16,3000,DPWAfile);
[fAr,~]=DPWAcpu(s,18,3000,DPWAfile);
[fCr,~]=DPWAcpu(s,24,3000,DPWAfile);
[fI,~]=DPWAcpu(s,53,3000,DPWAfile);

%% load molecule
fid=fopen(MoleculeFile); % the xyz file of the molecule
xyzFile=textscan(fid,'%s %f %f %f','headerlines', 2); % Skip the first two lines of .xyz file
fclose(fid);

AtomClass=xyzFile{1,1}; % the kind of different atoms
[AtomNum,~]=size(AtomClass); % the number of atoms
XYZpos=[xyzFile{1,2},xyzFile{1,3},xyzFile{1,4}]; % corresponding coordinates xyz 

% scattering amplitude
ScatAmp=[fH;fHe;fC;fN;fO;fF;fS;fAr;fCr;fI];

% to point the scattering amplitude in ScatAmp of atom i 
AtomPointer=zeros(1,AtomNum);
for ii=1:AtomNum
    switch char(AtomClass(ii))
        case 'H'
            AtomPointer(ii)=1;
        case 'He'
            AtomPointer(ii)=2;
        case 'C'
            AtomPointer(ii)=3;
        case 'N'
            AtomPointer(ii)=4;
        case 'O'
            AtomPointer(ii)=5;
        case 'F'
            AtomPointer(ii)=6;
        case 'S'
            AtomPointer(ii)=7;
        case 'Ar'
            AtomPointer(ii)=8;
        case 'Cr'
            AtomPointer(ii)=9;
        case 'I'
            AtomPointer(ii)=10;
        otherwise
            error('Such atom has not been added')
    end
end

%% calculating simulative IA and IM

SimIA=s.*0; % initialization IA=0
SimIM=SimIA; % initialization IM=0

for ii=1:AtomNum
    for jj=1:AtomNum
        if jj==ii
            SimIA=SimIA+abs(ScatAmp(AtomPointer(ii),:)).^2;
        else
            r=sqrt(sum((XYZpos(ii,:)-XYZpos(jj,:)).^2));
            SimIM=SimIM+ScatAmp(AtomPointer(ii),:).*conj(ScatAmp(AtomPointer(jj),:)).*sin(s.*r)./(s.*r);
        end 
    end
end
SimIM=real(SimIM);

SimIA=SimIA.*1e14;
SimIM=SimIM.*1e14;

end