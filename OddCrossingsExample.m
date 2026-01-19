% OddCrossingsExample.m
%

% S. Weiss, 30/3/2025

A = zeros(1,1,3);
A(1,1,:) = [1 0 1]/sqrt(2); 
B = zeros(1,1,1); B(1,1,1) = 1;
Q = PUProcrustes(A,B,1024,0,20);  
close all;

q = squeeze(Q);
stem(-256:255,q,'b');
hold on;
plot(0:2,squeeze(A),'r*');
plot((-10:-1),zeros(10,1),'r*'); plot(3:15,zeros(13,1),'r*');
axis([-10 15 -.25 .75]);
grid on;
xlabel('time index $n$','interpreter','latex');
ylabel('impulse response','interpreter','latex');  
legend({'$\hat{q}_{\ast}[n]$','$a[n]$'},'interpreter','latex','location','NorthEast');

set(gcf,'OuterPosition',[230 250 570 300]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc OddZeroCrossingsFig2.eps     




