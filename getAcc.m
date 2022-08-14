% Finding acceleration

function ax = getAcc(H,L,d,T,z,x,t)
a = H/2;
g = 9.81;
k = 2*pi/L;
w = 2*pi/T;

ax = (a*g*k*cosh(k*(d+z))/cosh(k*d))*sin(k*x - w*t);
end