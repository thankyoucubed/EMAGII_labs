% ECE332 lab1 part 2

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


theta1 = atan2(Y,(X+.005));           %theta sets the direction of phi from
theta2 = atan2(Y,(X-.005));           %cylindrical to cartesian coordinates

%Establishing the x and y components from both wires to be used in
%superposition for total B. Note r is replaced with it's cartesian
%equivalent of sqrt(x^2 + y^2) and reference of wires (the +/- .005 is also
%taken into account wrt the origin.

B1x = -(mu*I./(2*pi*sqrt((X+.005).^2 + Y.^2))).*sin(theta1);
B2x = -(mu*I./(2*pi*sqrt((X-.005).^2 + Y.^2))).*sin(theta2);

B1y = (mu*I./(2*pi*sqrt((X+.005).^2 + Y.^2))).*cos(theta1);
B2y = (mu*I./(2*pi*sqrt((X-.005).^2 + Y.^2))).*cos(theta2);

Bx = B1x + B2x;
By = B1y + B2y;

index1 = 1 : 20 : N;
index2 = index1;
          
p1 = 1.*X(index1, index2); 
p2 = 1.*Y(index1, index2);

Bscale = sqrt(Bx(index1,index2).^2 + By(index1,index2).^2);     
scale = 0.00001;
Bscale(Bscale < scale)= scale; 

p3 = 1e3.*Bx(index1, index2)./Bscale;
p4 = 1e3.*By(index1, index2)./Bscale;

figure(2)
f1=quiver(p1,p2,p3,p4,'autoscalefactor',0.5);
set(f1,'color',[0 0 1],'linewidth',1.2)
axis tight
title('B field of 2 Parallel Currents in Inf. Wire')
xlabel('x [m]')
ylabel('y [m]')

% Same deal but for opp. currents, chose wire 2 to be I = -10A (-z direction)

theta1 = atan2(Y,(X+.005));
theta2 = atan2(Y,(X-.005));

B1x = -(mu*I./(2*pi*sqrt((X+.005).^2 + Y.^2))).*sin(theta1);
B2x = -(mu*-I./(2*pi*sqrt((X-.005).^2 + Y.^2))).*sin(theta2);

B1y = (mu*I./(2*pi*sqrt((X+.005).^2 + Y.^2))).*cos(theta1);
B2y = (mu*-I./(2*pi*sqrt((X-.005).^2 + Y.^2))).*cos(theta2);

Bx = B1x + B2x;
By = B1y + B2y;

index1 = 1 : 20 : N;
index2 = index1;
          
p1 = 1.*X(index1, index2); 
p2 = 1.*Y(index1, index2);

Bscale = sqrt(Bx(index1,index2).^2 + By(index1,index2).^2);     
scale = 0.00001;
Bscale(Bscale < scale)= scale; 

p3 = 1e3.*Bx(index1, index2)./Bscale;
p4 = 1e3.*By(index1, index2)./Bscale;

figure(3)
f1=quiver(p1,p2,p3,p4,'autoscalefactor',0.5);
set(f1,'color',[0 0 1],'linewidth',1.2)
axis tight
title('B field of 2 Anti-Parallel Currents in Inf. Wire')
xlabel('x [m]')
ylabel('y [m]')

% Comments For Part 2
% The two parallel wires are expected to roughly cancel out (1/r distancing
% is seen) in the middle and look like they "couple" on the outside to make 
% one big flowing B vector field pointing in the CCW direction. The two 
% antiparallel wires are expected to look like reflections about the y axis
% as the the -I wire will now "couple" with the the +I wire in the center.

