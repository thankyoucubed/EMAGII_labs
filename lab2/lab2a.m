% Regan Garner
% Nicholas Porter
% Grace Semerjian

% ECE 332 Lab 2: Waveguides 

% make the following plots of Ez:

% (i) a 1-d plot of |Ez| vs. x at y=b/2

syms x z;

epsilon = 8.854*(10^-12); % electric permittivity / dielectric constant 
mu = 4*pi*(10^7); % magnetic permeability/ magnetic constant 

a = 1; % arbitrary width (x direction)
b = a/2; % height (y direction)
m = 1; % m and n are integers that indicate the mode of propagation of the wave.
n = 1; %
f = 20
omega = 2*pi*f;
kx = (m*pi)/a;
ky = (n*pi)/b;
kc = sqrt(kx^2+ky^2);
k = omega*sqrt(epsilon*mu); 
lambda = 3*10^8/f;

%condition for a is that
% a > (m*lambda)/2

% condition for b
% b > m*lambda

% cutoff freq
up0 = 1/sqrt(mu*epsilon);
fmn = up0/2*sqrt((m/a)^2+(n/b)^2)

beta = sqrt(k^2-kc^2); % the imaginary part of the propagation constant, the phase
% constant. Because we replace gamma with beta, we assume that the wave is unattenuated.

y = b/2; % fixed value for part (i)
z = 0;

Eo = 1; % arbitrary strength of the electric field 

% MATLAB is poopid, so we have to redefine the previously symbolic 
% variables as fixed constants in order to use fplot later 
E = Eo*sin((m*pi*x)/a)*sin((n*pi*y)/b)*exp(-1i*beta*z)

figure(1)
fplot(E)
title('1-d plot of |E| at y=b/2')
xlabel('x (m)');ylabel('|E| (V/m)')

% (ii) a 2-d contour plot of |Ez| in the xy plane
syms y
z = 1;
count=2;
for f=1:5:51
    omega = 2*pi*f;
    k = omega*sqrt(epsilon*mu);
    beta = sqrt(k^2-kc^2);
    [x,y] = meshgrid(0:.1:a,0:.05:b);
    Ex = ((-1i*beta)/kc^2)*(m*pi/a)*Eo*cos(m*pi*x/a)*sin(n*pi*y/b)*exp(-1i*beta*z);
    Ey = ((-1i*beta)/kc^2)*(n*pi/a)*Eo*sin(m*pi*x/a)*cos(n*pi*y/b)*exp(-1i*beta*z);
    figure(count)
    quiver(x,y,Ex,Ey);
    xlim([0 a]);ylim([0 b]);
    title('2-d contour plot of |E| in x-y plane for z=' + string(f))
    ylabel('y (m)');xlabel('x (m)');
    count = count+1;
end

% (iii) a 3-d surface plot of |Ez| in the xy plane
figure(count)
surf(Ex,Ey,z)
count = count +1;
