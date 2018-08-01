x = [0:0.01:1];
fx = -1 .* x.* x + 2 .* x;
figure(1);
plot (x, fx);

theta1 = 0.1
theta2 = 0.4
theta3 = 0.6
theta4 = 0.1

sigmaA = 0.1
sigmaB = 0.1
sigmaC = 0.2

kA =  0.1
kB =  6
kC =  0.9

rhoA = (kA./2) .* (1+tanh(fx-theta1)/sigmaA);
rhoB = (kB./2) .* (1+tanh((theta2-fx)/sigmaB)) .* (kB./2) .* (1+tanh((fx-theta3)./sigmaB));
rhoC = (kC./2) .* (1+tanh((theta4-fx)./sigmaC));

figure(2); clf;
hold on;
plot (x, rhoA, 'r')
plot (x, rhoB, 'b')
plot (x, rhoC, 'g')
