load ../logs/factorexpression.h5

figure (1);
clf;
plot3 (x,y,rhoA,'ro');
hold on;
plot3 (x,y,rhoB,'bo');
plot3 (x,y,rhoC,'go');
view([4.7,0]);
legend('\rho_A','\rho_B','\rho_C');
title('Signalling molecule expression')

figure (2);
clf;
plot3 (x,y,grad_rhoA_x,'ro');
hold on;
plot3 (x,y,grad_rhoB_x,'bo');
plot3 (x,y,grad_rhoC_x,'go');
view([4.7,0]);
legend('\rho_A','\rho_B','\rho_C');
title('Signalling molecule expression gradients')