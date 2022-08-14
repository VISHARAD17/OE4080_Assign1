% Finding force

function force = getForce(Cd,Cm,D,H,L,z,d,x,T,t)
rho = 1025;
vel = getVelocity(H,L,d,T,z,x,t);
acc = getAcc(H,L,d,T,z,x,t);
V = pi*D^2/4;
Fi = rho*Cm*V*acc;
Fd = 0.5*rho*Cd*D*vel*abs(vel);

% Force per unit length can be calculated by
force = Fi + Fd;
end