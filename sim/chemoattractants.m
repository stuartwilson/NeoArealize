%
% A script to find good parameters for sigmaA-C, kA-C and theta1-4
%

load ../logs/factorexpression.h5

% Use a loaded version of fgf
idxs = find (y(find (y>-0.2)) < 0.2);
x_ = x(idxs);
fx = fgf(idxs);
figure(10);
plot (x_, fx, 'o');

sigmaA = 0.2
sigmaB = 0.2
sigmaC = 0.2

kA =  0.7
kB =  0.9
kC =  0.48

theta1 = 0.9
theta2 = 0.5
theta3 = 0.39
theta4 = 0.08

_rhoA = (kA./2) .* (1+tanh((fx-theta1)./sigmaA));
_rhoB = (kB./2) .* (1+tanh((theta2-fx)./sigmaB)) .* (kB./2) .* (1+tanh((fx-theta3)./sigmaB));
_rhoC = (kC./2) .* (1+tanh((theta4-fx)./sigmaC));

figure(11); clf;
hold on;
plot (x_, _rhoA, 'ro')
plot (x_, _rhoB, 'bo')
plot (x_, _rhoC, 'go')
%plot (x_, rhoA(idxs), 'mo');