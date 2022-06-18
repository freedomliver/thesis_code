function [nonspinflip,spinflip] = DPWAcpu(q,atom,Ek,DPWAfile)
% DPWA (Dirac partial-wave analysis; calculation of atomic elastic scattering of electrons)
% using cpu
% 2020/8/31 qifengfeng
% 
% q: the the momentum transfer related to scattering angular;
%      q=2ksin(theta/2); an array
%      (angstrom^-1)
% atom: atomic number; the kind of atom; int
% Ek: the kinetic energy of incident electrons (keV); single value
% DPWAfile: the path of the file DPWA
%
% nonspinflip: f(q) the direct scattering amplitude (cm); an array
% spinflip: g(q) spin-flip scattering amplitude (cm); an array

% international unit
me=9.10938356e-31;  % mass of electron
qe=1.602176634e-19;  % elementary charge
Eke=Ek*1e3*qe;  % kinetic energy of electron
hbar=6.626075540e-34/(2*pi);  % reduced Planck constant
c=2.99792458e8; % speed of light
k=sqrt(Eke*(Eke+2*me*c^2))/(c*hbar);  % wave vector of electron

q=q.*1e10; % (m^-1)
theta=2.*asin(q./2./k);
input=cos(theta);

% load DPWA 
DPWA_file=strcat('dpwa',num2str(atom),'_',num2str(Ek),'.mat');
path=strcat(DPWAfile,DPWA_file); % DPWA data from ELSEPA saved in the file
if exist(path,'file')==0
    error('Such DPWA has not been calculated from ELSEPA or the input is not valid or the path of DPWA file is not correct.')
end
tmp=load(path);
names = fieldnames(tmp);
DPWA=tmp.(names{1});

% the direct scattering amplitude; f(theta)
l=0; % orbital quantum number l
spinup=DPWA(1,2); % phase shift of s=1/2
spindown=DPWA(1,3); % phase shift of s=-1/2
tmp=size(q);
P1=ones(tmp); % legendre function P0(x)=1;
nonspinflip=((l+1)*(exp(2*1i*spinup)-1)+l*(exp(2*1i*spindown)-1)).*P1;  % f(theta)

l=1; 
spinup=DPWA(2,2); % phase shift of s=1/2
spindown=DPWA(2,3); % phase shift of s=-1/2
P2=input; % legendre function P1(x)=x;
nonspinflip=nonspinflip+((l+1)*(exp(2*1i*spinup)-1)+l*(exp(2*1i*spindown)-1)).*P2;

tmp=size(DPWA);
l_num=tmp(1);  % the number of orbital quantum number l in DPWA
for ll=3:l_num % from l=2
    l=DPWA(ll,1);
    spinup=DPWA(ll,2); % phase shift of s=1/2
    spindown=DPWA(ll,3); % phase shift of s=-1/2
    P3=((2*l-1).*input.*P2-(l-1).*P1)./l; % legendre recurrence for m=0
    P1=P2;
    P2=P3;
    nonspinflip=nonspinflip+((l+1)*(exp(2*1i*spinup)-1)+l*(exp(2*1i*spindown)-1)).*P3;
end

% spin-flip scattering amplitude; g(theta)
% l=1; % orbital quantum number l
spinup=DPWA(2,2); % phase shift of s=1/2
spindown=DPWA(2,3); % phase shift of s=-1/2
Q1=-(1-input.^2).^(1/2); % legendre function P11(x)=-(1-x^2)^(1/2)
spinflip=(exp(2*1i*spindown)-exp(2*1i*spinup)).*Q1; % g(theta)

% l=2; 
spinup=DPWA(3,2); % phase shift of s=1/2
spindown=DPWA(3,3); % phase shift of s=-1/2
Q2=-3.*input.*(1-input.^2).^(1/2); % legendre function P21(x)=-3x(1-x^2)^(1/2);
spinflip=spinflip+(exp(2*1i*spindown)-exp(2*1i*spinup)).*Q2;

for ll=4:l_num % from l=3
    l=DPWA(ll,1);
    spinup=DPWA(ll,2); % phase shift of s=1/2
    spindown=DPWA(ll,3); % phase shift of s=-1/2
    Q3=((2*l-1).*input.*Q2-l.*Q1)./(l-1); % legendre recurrence for m=1
    Q1=Q2;
    Q2=Q3;
    spinflip=spinflip+(exp(2*1i*spindown)-exp(2*1i*spinup)).*Q3;
end

nonspinflip=nonspinflip./(2*1i*k).*1e2; % (cm)
spinflip=spinflip./(2*1i*k).*1e2; % (cm)

end