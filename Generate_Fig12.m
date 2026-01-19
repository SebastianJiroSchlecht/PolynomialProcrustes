% Generate Figure 12 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., submitted Mar. 2025
%
% example for frequency domain of solution

clear all; close all;

if exist('Generate_Fig12.mat')~=2,
% matrix definition
S = zeros(2,2,3);
S(1,1,:) = [0 1 0];
S(2,2,:) = [.5 .5 .5];

U = zeros(2,2,2);
U(:,:,1) = [1 1; 0 0]/sqrt(2);
U(:,:,2) = [0 0; 1 -1]/sqrt(2);
V = zeros(2,2,1);
V(:,:,1) = [1 -1; 1 1]/sqrt(2);
A3 = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

A = zeros(2,2,4);
A(1,1,:) = [-1 1 -1 0]/4;
A(1,2,:) = [1 3 1 0]/4;
A(2,1,:) = [0 1 3 1]/4;
A(2,2,:) = [0 -1 1 -1]/4;

% DFT
Af = fft(A,1024,3);

% binwise solution (see Sebastian)
Qfb = zeros(2,2,1024);
for i = 1:1024,
   [u,~,v] = svd(Af(:,:,i));
   Qfb(:,:,i) = u*v';
end;    

% analytic solution
N = 16;                                      % this is the order of the allpass filter --- feel free to vary
A2 = zeros(2,2,size(A3,3)+N+1);
A2(:,:,N+2:end) = A3;
B = zeros(2,2,1); B(:,:,1) = eye(2);
Qhat = PUProcrustes(A2,B,2^12,0,N);  
ProcrustesMetrics(A3,S,U,V,Qhat);  

Nfft = size(Qhat,3);
Qhat2 = zeros(2,2,4*Nfft);
Qhat2(:,:,1:Nfft) = Qhat;
Nfft = 4*Nfft;
Qf = fft(circshift(Qhat2,[0 0 -N-1]),Nfft,3);

% display of results
disp('A --- solid lines, with blue real and red imaginary');
disp('analytic Q --- dashed lines with asterisks');
disp('binwise Q --- dash-dotted lines with circles'); 
   save Generate_Fig12.mat
else
   load Generate_Fig12.mat
end;
   
figure(8); clf;
ri = 2; ci = 2; FS = 10;
Seg1 = (1:342); Seg2 = (343:683); Seg3 = (684:1024);
f = (0:1023)/1024;

% for legend only
h = plot([-1 -1],[-1 -1],'-','linewidth',3); set(h(1),'color',[1 1 1]*.75);
hold;
plot([-1 -1],[-1 -1],'k--');
plot([-1 -1],[-1 -1],'k-');

%------------------------------------------------------
% plot main curves
%------------------------------------------------------
h = plot((0:1023)/1024,real(squeeze(Af(ri,ci,:))),'-','linewidth',3); 
set(h(1),'color',[.8 0.8 1]);
plot((0:Nfft-1)/Nfft,real(squeeze(Qf(ri,ci,:))),'b-'); 
plot(f(Seg1),real(squeeze(Qfb(ri,ci,Seg1))),'b--');
h1 = plot(f(Seg1(end)),real(squeeze(Qfb(ri,ci,Seg1(end)))),'bo','LineWidth',1.5); 
set(h1, 'markerfacecolor','r'); 
plot(f(Seg2),real(squeeze(Qfb(ri,ci,Seg2))),'b--'); 
h1 = plot(f(Seg2(1)),real(squeeze(Qfb(ri,ci,Seg2(1)))),'bo','LineWidth',1.5); 
set(h1, 'markerfacecolor','b'); 
h1= plot(f(Seg2(end)),real(squeeze(Qfb(ri,ci,Seg2(end)))),'bo','LineWidth',1.5); 
%set(h1, 'markerfacecolor','b'); 
h1 = plot(f(Seg3(1)),real(squeeze(Qfb(ri,ci,Seg3(1)))),'ro','LineWidth',1.5); 
set(h1, 'markerfacecolor','b'); 
plot(f(Seg3),real(squeeze(Qfb(ri,ci,Seg3))),'b--'); 
h = plot((0:1023)/1024,imag(squeeze(Af(ri,ci,:))),'-','linewidth',3); 
set(h(1),'color',[1 0.6 0.6]);
plot((0:Nfft-1)/Nfft,imag(squeeze(Qf(ri,ci,:))),'r-');
plot(f(Seg1),imag(squeeze(Qfb(ri,ci,Seg1))),'r--');
plot(f(Seg2),imag(squeeze(Qfb(ri,ci,Seg2))),'r--');
plot(f(Seg2(1)),imag(squeeze(Qfb(ri,ci,Seg2(1)))),'ro','LineWidth',1.5);
h1 = plot(f(Seg2(end)),imag(squeeze(Qfb(ri,ci,Seg2(end)))),'ro','LineWidth',1.5);
set(h1, 'markerfacecolor','r'); 
plot(f(Seg3),imag(squeeze(Qfb(ri,ci,Seg3))),'r--');
axis([0 1 -1.1 1.1]);
% boxes to indicate inserts
b = [23.5/36 24.5/36  -1 .2];
boxx1 = [b(1) b(2) b(2) b(1) b(1)]; boxy1 = [b(4) b(4) b(3) b(3) b(4)];
plot(boxx1,boxy1,'k');
b = [11.5/36 12.5/36 -.8 1];
boxx1 = [b(1) b(2) b(2) b(1) b(1)]; boxy1 = [b(4) b(4) b(3) b(3) b(4)];
plot(boxx1,boxy1,'k');

set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/6:1),'XTickLabel',{'$0$','$\pi/3$','$2\pi/3$','$\pi$','$4\pi/3$',...
      '$5\pi/3$','$2\pi$'});
%    'YTick',(0:1:2),'YTickLabel',...
%     {'$0$','$1$','$2$'});
legend({'$A_{22}(\mathrm{e}^{\mathrm{j}\Omega})$',...
       '$\hat{Q}_{22}^{(N \rightarrow \infty)}(\Omega)$',...
       '$\hat{Q}_{\ast,22}(\mathrm{e}^{\mathrm{j}\Omega})$'},...
        'interpreter','latex','location','SouthEast');
grid on;        
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
ylabel('$\Re\{\cdot\}$, $\Im\{\cdot\}$',...
	'interpreter','latex','fontsize',FS);
        
%------------------------------------------------------
% plot upper right insert
%------------------------------------------------------
axes('position',[.71 .65 .18 .25]);
plot((0:Nfft-1)/Nfft,real(squeeze(Qf(ri,ci,:))),'b-'); 
hold on;
plot(f(Seg1),real(squeeze(Qfb(ri,ci,Seg1))),'b--');
h1 = plot(f(Seg1(end)),real(squeeze(Qfb(ri,ci,Seg1(end)))),'bo','LineWidth',1.5); 
set(h1, 'markerfacecolor','r'); 
plot(f(Seg2),real(squeeze(Qfb(ri,ci,Seg2))),'b--'); 
h1 = plot(f(Seg2(1)),real(squeeze(Qfb(ri,ci,Seg2(1)))),'bo','LineWidth',1.5); 
set(h1, 'markerfacecolor','b'); 
h1= plot(2/3,real(squeeze(Qfb(ri,ci,Seg2(end)))),'bo','LineWidth',1.5); 
%set(h1, 'markerfacecolor','b'); 
h1 = plot(2/3,real(squeeze(Qfb(ri,ci,Seg3(1)))),'ro','LineWidth',1.5); 
set(h1, 'markerfacecolor','b'); 
plot(f(Seg3),real(squeeze(Qfb(ri,ci,Seg3))),'b--'); 
%h = plot((0:1023)/1024,imag(squeeze(Af(ri,ci,:))),'-','linewidth',3); 
%set(h(1),'color',[1 0.6 0.6]);
plot((0:Nfft-1)/Nfft,imag(squeeze(Qf(ri,ci,:))),'r--');
plot(f(Seg1),imag(squeeze(Qfb(ri,ci,Seg1))),'r--');
plot(f(Seg2),imag(squeeze(Qfb(ri,ci,Seg2))),'r--');
plot(f(Seg2(1)),imag(squeeze(Qfb(ri,ci,Seg2(1)))),'ro','LineWidth',1.5);
h1 = plot(2/3,imag(squeeze(Qfb(ri,ci,Seg2(end)))),'ro','LineWidth',1.5);
set(h1, 'markerfacecolor','r'); 
plot(f(Seg3),imag(squeeze(Qfb(ri,ci,Seg3))),'r--');
plot((0:Nfft-1)/Nfft,imag(squeeze(Qf(ri,ci,:))),'r-');
axis([23.5/36 24.5/36 -1 .2]);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(23.5:.5:24.5)/36,'XTickLabel',{'$\frac{47\pi}{36}$',...
     '$\frac{4\pi}{3}$','$\frac{49\pi}{36}$'});
grid on;

%------------------------------------------------------
% plot upper left insert
%------------------------------------------------------
axes('position',[.175 .65 .18 .25]);
plot((0:Nfft-1)/Nfft,real(squeeze(Qf(ri,ci,:))),'b-'); 
hold on;
plot(f(Seg1),real(squeeze(Qfb(ri,ci,Seg1))),'b--');
h1 = plot(1/3,real(squeeze(Qfb(ri,ci,Seg1(end)))),'bo','LineWidth',1.5); 
set(h1, 'markerfacecolor','r'); 
plot(f(Seg2),real(squeeze(Qfb(ri,ci,Seg2))),'b--'); 
h1 = plot(1/3,real(squeeze(Qfb(ri,ci,Seg2(1)))),'bo','LineWidth',1.5); 
set(h1, 'markerfacecolor','b'); 
h1= plot(f(Seg2(end)),real(squeeze(Qfb(ri,ci,Seg2(end)))),'bo','LineWidth',1.5); 
%set(h1, 'markerfacecolor','b'); 
h1 = plot(f(Seg3(1)),real(squeeze(Qfb(ri,ci,Seg3(1)))),'ro','LineWidth',1.5); 
set(h1, 'markerfacecolor','b'); 
plot(f(Seg3),real(squeeze(Qfb(ri,ci,Seg3))),'b--'); 
%h = plot((0:1023)/1024,imag(squeeze(Af(ri,ci,:))),'-','linewidth',3); 
%set(h(1),'color',[1 0.6 0.6]);
plot((0:Nfft-1)/Nfft,imag(squeeze(Qf(ri,ci,:))),'r--');
plot(f(Seg1),imag(squeeze(Qfb(ri,ci,Seg1))),'r--');
plot(f(Seg2),imag(squeeze(Qfb(ri,ci,Seg2))),'r--');
plot(1/3,imag(squeeze(Qfb(ri,ci,Seg2(1)))),'ro','LineWidth',1.5);
h1 = plot(1/3,imag(squeeze(Qfb(ri,ci,Seg2(end)))),'ro','LineWidth',1.5);
set(h1, 'markerfacecolor','r'); 
plot(f(Seg3),imag(squeeze(Qfb(ri,ci,Seg3))),'r--');
plot((0:Nfft-1)/Nfft,imag(squeeze(Qf(ri,ci,:))),'r-');
axis([11.5/36 12.5/36 -.8 1]);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(11.5:.5:12.5)/36,'XTickLabel',{'$\frac{23\pi}{36}$',...
     '$\frac{2\pi}{3}$','$\frac{25\pi}{36}$'});
     grid on;
     
set(gcf,'OuterPosition',[230 250 570 350]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig12.eps     

