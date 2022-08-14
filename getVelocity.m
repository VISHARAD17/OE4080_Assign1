% Finding velocity

function u = get_velocity(H,L,d,T,z,x,t)
w = 2*pi/T;
k = 2*pi/L;
u = pi*H*cosh(k*(d+z))*cos(k*x - w*t)/T/sinh(k*d);
end