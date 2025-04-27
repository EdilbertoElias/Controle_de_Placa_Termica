function [G0,T1,L] = parametrosFOPTD(h,y,L,Ts)

R = zeros(3);
f = zeros(3,1);
N = length(y);

deltat = Ts; inicio = L/deltat;
for k = inicio:N
    phi = [h*(deltat*k) -h -y(k)]';
    R = R + phi*phi';
    A = sum(y(1:k))*deltat;
    f = f + phi*A;
end

theta = R\f;
G0 = theta(1);
L = theta(2)/G0;
T1 = theta(3);