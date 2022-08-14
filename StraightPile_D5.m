clear all;
close all;
clc;

% 0 for extrapolated airy wave theory
% 1 for wheeler's stretching
% 2 for both
method = 2;

%%% Given Data
H = 18;
a = H/2;
d = 85;
T = 14;
D = 5;
L = getWavelength(d,T);
k = 2*pi/L;
w = 2*pi/T;
x = 0;
Cd = 0.7;
Cm = 2;
t_max = 20;
dt = linspace(0,t_max,1000);
Force_a = zeros(1,numel(dt));
Force_w = Force_a;
numel_sec = 1000;
fi_airy = zeros(numel_sec, numel(dt));
fi_wheeler = fi_airy;
Z = fi_airy; 

for i = 1:numel(dt)
    t = dt(i);
    eta = a*cos(k*x - w*t);
    z = linspace(eta,-d,1000);
    Z(:,i) = z;
    f = zeros(1,numel(z));
    
    %%%%% extrapolated airy wave theory
    if method == 0
        for j = 1:numel(z)
            if eta > 0 && z(j) > 0
                f(j) = geForce(Cd,Cm,D,H,L,0,d,x,T,t);
            else
                f(j) = getForce(Cd,Cm,D,H,L,z(j),d,x,T,t);
            end       
        end
        fi_airy(:,i) = f;
        Force_a(i) = trapz(z,f);

    %%%%% wheeler's stretching
    elseif method == 1   
        z_p = (z-eta)*(d/(d+eta));
        for j= 1:numel(z_p)
            f(j) = getForce(Cd,Cm,D,H,L,z_p(j),d,x,T,t);
        end
        fi_wheeler(:,i) = f;
    
        Force_w(i) = trapz(z_p,f); 

    elseif method == 2  % plot both
        %%% airy
        for j = 1:numel(z)
            if eta > 0 && z(j) > 0
                f(j) = getForce(Cd,Cm,D,H,L,0,d,x,T,t);
            else
                f(j) = getForce(Cd,Cm,D,H,L,z(j),d,x,T,t);
            end       
        end
        fi_airy(:,i) = f;
        Force_a(i) = trapz(z,f);
        
        %%% Wheeler
        z_p = (z-eta)*(d/(d+eta));
        for j= 1:numel(z_p)
            f(j) = getForce(Cd,Cm,D,H,L,z_p(j),d,x,T,t);
        end
        fi_wheeler(:,i) = f;
        Force_w(i) = trapz(z_p,f); 
    end
end

%%% Graph
if method == 0
    figure;
    hold on;
    plot(dt, Force_a, '.','LineWidth',2);
    title('Waveloads estimated using Extapolated Airy Wave Theory');
    xlabel('t (sec)')
    ylabel('Force (N)')
    grid on;
    hold off;
elseif method == 1
    figure;
    hold on;
    plot(dt, Force_w, '.','LineWidth',2);
    title('Waveloads estimated using Wheeler''s stretching');
    xlabel('t (sec)')
    ylabel('Force (N)')
    grid on;
    hold off;
elseif method == 2
    figure;
    hold on;
    plot(dt, Force_a, '.','LineWidth',2);
    plot(dt, Force_w, '.','LineWidth',2);
    legend('Extrapolated Airy Wave Theory','Wheeler''s Stretching')
    title('Plot for D = 5');
    xlabel('t (sec)')
    ylabel('Force (N)')
    grid on;
    hold off;
end
