%% The Title
%-------------------------------------------------------------------------%
%----------------Lab 1: Faraday's Law and Biot-Savart's Law---------------%
%-------------------------ver. 1, 03/28/2020------------------------------%
%--------CE332 Lab, Spring 2020: Special remote education edition---------%


%% Part 1: Faraday's Law

% A 10cm x 10cm square loop of wire with internal resistance 1 ohm has a 
% small gap in one edge closest to a coplanar long straight current-
% carrying wire starting from a distance 5cm away. The current is time-
% varying, described by a cosine with peak magnitude 5 A and frequency 
% f. Using Matlab, find the flux through the loop as a symbolic
% function of time for symbolic frequency f, then plot the induced EMF in 
% the loop if for f stepped by octave from 10 Hz to 160 Hz. Assume the 
% current carrying wire is the z axis, pointed in the +z direction, and 
% that all vertices of the loop of wire are in the y-z plane. Plot EMF over
% the longest period under consideration. In a separate section, comment on
% the direction of the flux and the polarity of the induced EMF: if there 
% were a 4 ohm resistance bridging the gap in the loop would the induced 
% current flow clockwise or counterclockwise from the perspective of a 
% viewer looking in the direction (-x)?

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
f = f_low:2*flow:f_hi;
t = linspace(0,1/f_low,length(f));

I(f,t) = 5*cos(2*pi.*f.*t);
B(y,f,t) = mu_0.*I(f,t)/(2*pi.*y);
flux(f,t) = int(int(B(y,f,t),z,0,.1),y,.05,.15);

fprintf('The total flux through the loop in Webers is: ')
disp(vpa(flux(f,t),4))

% Part 2: Calculate and display the Vemf across the small gap in the loop
% as a symbolic function emf(f,t).

emf(f,t) = -N.*diff(flux(f,t),t);

fprintf('The induced EMF in volts is: ')
disp(vpa(emf(f,t),4))

% Plot emf(t) for frequencies stepped by octave from 10 Hz to 160 Hz

plot(f,emf(f,t))

%% Comments For Part 1




%% Part 2: Ampere's Law and Vector Field Plotting

% Using a quiver plot, show the magnetic field B between two infinitely
% long parallel conductors 1cm apart each carrying 10 A, first as if they
% both carry current in the same direction (+z) and again when one is
% opposing (I in -z direction). 

% Ignore the effects of mutual inductive coupling and use superposition to 
% obtain the total field vectors at all points. 
 
clc
close all
clear all

N = 401;                                % 401 points per axis
I = 10;
dist = 0.01;
mu = 4*pi*1e-7;

x = linspace(-0.01,0.01,N);               % N points for X
y = linspace(-0.01,0.01,N);               % N points for y
[X,Y] = meshgrid(x,y);                  % create 2D coordinate space

wireX = [-dist/2 dist/2];            % X-coordinates of wires 1 and 2
wireY = [0 0];                       % Y-coordinates of wires 1 and 2


% YOUR CODE GOES HERE. PROVIDE A 401x401 MATRIX FOR EACH OF SETS, X AND Y 
% VECTORS OF THE B-FIELD, CALLED 'Bx' AND 'By'. THE BELOW CODE WILL
% SCALE IT TO MAKE FOR A NICE QUIVER PLOT. COMMENT ON WHAT THESE SNIPPETS
% ARE DOING AND CONSIDER WHY THAT MIGHT BE NECESSARY OR USEFUL IN THE
% SECTION BELOW.


index1 = 1 : 20 : N;
index2 = index1;
          
p1 = 1.*X(index1, index2); 
p2 = 1.*Y(index1, index2);

Bscale = sqrt(Bx(index1,index2).^2 + By(index1,index2).^2);     
scale = 0.00001;
Bscale(Bscale < scale)= scale; 

p3 = 1e3.*Bx(index1, index2)./Bscale;
p4 = 1e3.*By(index1, index2)./Bscale;


figure
f1=quiver(p1,p2,p3,p4,'autoscalefactor',0.5);
set(f1,'color',[0 0 1],'linewidth',1.2)
axis tight


%end

%% Comments For Part 2



