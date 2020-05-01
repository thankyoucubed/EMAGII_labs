% lab 3

c = 2.99792458*10^8;% speed of light
I0 = 1;             % current amplitude
f = 2.4*10^6;           % frequency
lambda = c/f;       % wavelength
k = (2*pi)/lambda;  % wave #
L = 1.25*lambda;    % length of antenna
N = 100;            % # of data points
z = linspace(-L/2,L/2,N);

% piecewise equation for current phasor I~(z)
%{
for i = 1:N
    if z(i) < 0
        I(i) = I0*sin(k*(L/2+z(i)));
    elseif abs(z(i)) > L/2
        disp('error: z out of bounds')
    else
        I(i) = I0*sin(k*(L/2-z(i)));
    end
end
%}
I = I0.*sin(k.*(L/2-abs(z)));

plot(z,I)