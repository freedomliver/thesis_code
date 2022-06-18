function [output] = FourierHankel(pattern,n,rnum)
% inverse Fourier transformation followed by inverse Hankel transformation
%
% pattern: pattern to be processed
% n: pattern size(length)
% rnum: radius of output data
%
% output: output data


% ifft
fft_pattern=pattern;
for i=1:n
    fft_pattern(:,i)=ifft(fftshift(pattern(:,i)));
    fft_pattern(:,i)=fftshift(fft_pattern(:,i));
end

% hankel inversion
N=n/2;
X=40;
[~,~,~,I,K,R]=dht([],X,N,0); 
data=zeros(n,N);
for j=1:n
    data(j,:)=idht(real(fft_pattern(j,N+1:n)'),I,K,R)';
end

% output data
output=[fliplr(data(n/2+1-rnum:n/2+1+rnum,2:rnum+1)),data(n/2+1-rnum:n/2+1+rnum,1:rnum+1)];

end