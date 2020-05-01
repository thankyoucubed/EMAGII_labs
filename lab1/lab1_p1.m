% ECE332 lab1 part 1

clc
close all
clear all

syms t f x y z;
mu_0 = pi*4e-7;

% Part 1)
% Calculate the Flux through the area of the loop. Round the expression to
% 4 significant figures (use the vpa() function for this).

N = 1;
f_low = 10;
f_hi = 160;

I(f,t) = 5*cos(2*pi*f*t);
B(y,f,t) = mu_0*I(f,t)/(2*pi*y);
flux(f,t) = int(int(B(y,f,t),z,0,.1),y,.05,.15);

fprintf('The total flux through the loop in Webers is: ')
disp(vpa(flux(f,t),4))

% Part 2: Calculate and display the Vemf across the small gap in the loop
% as a symbolic function emf(f,t).

emf(f,t) = -N*diff(flux(f,t),t);

fprintf('The induced EMF in volts is: ')
disp(vpa(emf(f,t),4))

% Plot emf(t) for frequencies stepped by octave from 10 Hz to 160 Hz
f = f_low;
t = linspace(0,1/f_low);

plot(t,emf(f,t))
hold on
for k=1:5
    f = 2*f
    plot(t,emf(f,t))
end
hold off
title('Induced V_e_m_f in the loop for f stepped by octave from 10Hz to 160Hz')
xlabel('time (s)');ylabel('V_e_m_f');