function p = LegendreZero(n,x)
% Legendre polynomials, m=0
%
% n: Legendre order
% x: cos(theta)
%
% p: output data

if (n==0)
    p=1;
else
    if (n==1)
        p=x;
    else
        p=((2.*n-1).*x.*Legendre(n-1,x)-(n-1).*Legendre(n-2,x))./n;
    end
end

end
