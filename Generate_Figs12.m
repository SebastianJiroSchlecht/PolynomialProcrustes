% Generate Figures 1 and 2 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., submitted Mar. 2025
%
% showcase the purpose of the ultimate Procrustes solution

clear all; close all;

S = zeros(2,2,3);
S(1,1,:) = [0 1 0];
S(2,2,:) = [.5 .5 .5];

s1 = zeros(1,1024); s1(1:2) = [1 0]; s1(end) = 0;
s2 = zeros(1,1024); s2(1:2) = [.5 .5]; s2(end) = .5;
U = zeros(2,2,2);
U(:,:,1) = [1 1; 0 0]/sqrt(2);
U(:,:,2) = [0 0; 1 -1]/sqrt(2);
V = zeros(2,2,1);
V(:,:,1) = [1 1; -1 1]/sqrt(2);
A = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

AAP = PolyMatConv(A,ParaHerm(A));
I = zeros(2,2,1); I(:,:,1)=eye(2);

[Q,D,Sfcorr2,OSflag,Nfft,e] = PUProcrustes(A,I,1024,0,32);

%-------------------------------------------------------------------------------
% Figure 1 --- example matrix and Procrustes solution
%-------------------------------------------------------------------------------
figure(1); clf;
FS = 12;
Anew = zeros(2,2,26);
Anew(:,:,11:14) = A;
Ahat3 = real(Q(:,:,503:528));
t = (-10:15); FS= FS-2;
subplot(221);
stem(t,squeeze(Anew(1,1,:)),'b'); hold on;
plot(t,squeeze(Ahat3(1,1,:)),'r*');
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
ylabel('$a_{1,m}[n], \; \hat{q}_{1,m}[n]$','interpreter','latex','fontsize',FS);
text(-9,0.75,'$m=1$','interpreter','latex','fontsize',FS);
subplot(222);
stem(t,squeeze(Anew(1,2,:)),'b'); hold on;
plot(t,squeeze(Ahat3(1,2,:)),'r*');
text(-9,0.75,'$m=2$','interpreter','latex','fontsize',FS);
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
subplot(223);
stem(t,squeeze(Anew(2,1,:)),'b'); hold on;
plot(t,squeeze(Ahat3(2,1,:)),'r*');
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
text(-9,0.75,'$m=1$','interpreter','latex','fontsize',FS);
ylabel('$a_{2,m}[n], \; \hat{q}_{2,m}[n]$','interpreter','latex','fontsize',FS);
xlabel('time index $n$','interpreter','latex','fontsize',FS);
subplot(224);
stem(t,squeeze(Anew(2,2,:)),'b'); hold on;
plot(t,squeeze(Ahat3(2,2,:)),'r*');
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
text(-9,0.75,'$m=2$','interpreter','latex','fontsize',FS);
xlabel('time index $n$','interpreter','latex','fontsize',FS);
set(gcf,'OuterPosition',[230 250 570 350]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig1.eps

%-------------------------------------------------------------------------------
% Figure 2 --- show AA^P and QQ^P
%-------------------------------------------------------------------------------
figure(2); clf;
FS = 12;
RA = PolyMatConv(A,ParaHerm(A));
RQ = real(PolyMatConv(Q,ParaHerm(Q)));
Ahat3 = RQ(:,:,1014:1034);
Anew = zeros(2,2,21);
Anew(:,:,8:14) = RA;
t = (-10:10); FS= FS-2;
subplot(221);
stem(t,squeeze(Anew(1,1,:)),'b'); hold on;
plot(t,squeeze(Ahat3(1,1,:)),'r*');
axis([-10.5 10.5 -.5 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
ylabel('$r_{A,1,m}[n], \; r_{Q,1,m}[n]$','interpreter','latex','fontsize',FS);
text(-9,0.75,'$m=1$','interpreter','latex','fontsize',FS);
subplot(222);
stem(t,squeeze(Anew(1,2,:)),'b'); hold on;
plot(t,squeeze(Ahat3(1,2,:)),'r*');
text(-9,0.75,'$m=2$','interpreter','latex','fontsize',FS);
axis([-10.5 10.5 -.5 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
subplot(223);
stem(t,squeeze(Anew(2,1,:)),'b'); hold on;
plot(t,squeeze(Ahat3(2,1,:)),'r*');
axis([-10.5 10.5 -.5 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
text(-9,0.75,'$m=1$','interpreter','latex','fontsize',FS);
ylabel('$r_{A,2,m}[n], \; r_{Q,2,m}[n]$','interpreter','latex','fontsize',FS);
xlabel('time index $n$','interpreter','latex','fontsize',FS);
subplot(224);
stem(t,squeeze(Anew(2,2,:)),'b'); hold on;
plot(t,squeeze(Ahat3(2,2,:)),'r*');
axis([-10.5 10.5 -.5 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
text(-9,0.75,'$m=2$','interpreter','latex','fontsize',FS);
xlabel('time index $n$','interpreter','latex','fontsize',FS);
set(gcf,'OuterPosition',[230 250 570 350]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig2.eps
