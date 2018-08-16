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

figure (3);
clf;
plot3 (x,y,fgf,'ro');
view([4.7,0]);
legend('fgf');
title('fgf expression')

figure (4);
clf;
plot3 (x,y,g_0_x,'ro');
view([4.7,0]);
zlim([-8 8]);
legend('g\_0_x');
title('g modifier')

figure (5);
clf;
plot3 (x,y,g_1_x,'yo');
view([4.7,0]);
zlim([-8 8]);
legend('g\_1_x');
title('g modifier')

figure (6);
clf;
plot3 (x,y,g_2_x,'go');
view([4.7,0]);
zlim([-8 8]);
legend('g\_2_x');
title('g modifier')

figure (7);
clf;
plot3 (x,y,g_3_x,'bo');
view([4.7,0]);
zlim([-8 8]);
legend('g\_3_x');
title('g modifier')

figure (8);
clf;
plot3 (x,y,g_4_x,'mo');
view([4.7,0]);
zlim([-8 8]);
legend('g\_4_x');
title('g modifier')
