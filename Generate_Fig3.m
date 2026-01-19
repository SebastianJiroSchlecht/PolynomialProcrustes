% Generate Figure 3 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., submitted Mar. 2025
%
% Example for an analytic SVD with a negative singular value

clear all; close all;

FS = 12;

% singular values
S = zeros(2,2,3);
S(1,1,:) = [0 1 0];
S(2,2,:) = [.5 .5 .5];

%------------------------------------------------------------------------------
%  Figure 3 --- singular values, Fourier domain
%------------------------------------------------------------------------------
s1 = zeros(1,1024); s1(1:2) = [1 0]; s1(end) = 0;
s2 = zeros(1,1024); s2(1:2) = [.5 .5]; s2(end) = .5;
plot((0:1023)/1024,real(fft(s1,1024)),'b-');
hold on;
plot((0:1023)/1024,real(fft(s2,1024)),'r--');
axis([0 1 -.5 1.5]); grid on;
ylabel('$\sigma_m(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/6:1),'XTickLabel',{'$0$','$\pi/3$','$2\pi/3$','$\pi$','$4\pi/3$','$5\pi/3$',...
      '$2\pi$'},...
    'YTick',(-1.5:.5:1.5),'YTickLabel',...
     {'$-1.5$','$-1$','$-.5$','$0$','$.5$','$1$','$1.5$'});
legend({'$m=1$','$m=2$'},'interpreter','latex','fontsize',FS-2,...
        'location','West');
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
set(gcf,'OuterPosition',[230 250 570 250]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig.eps



