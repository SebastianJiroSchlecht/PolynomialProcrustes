% Generate Figures 4, 5, 6, and 7 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., submitted Mar. 2025
%

clear all; close all;

FS = 12;

%------------------------------------------------------------------------------
%  systems and parameters
%------------------------------------------------------------------------------
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

%A2 = zeros(2,2,4);
%A2(:,:,1) = [-1 1; 0 0];
%A2(:,:,2) = [1 3;  1 -1];
%A2(:,:,3) = [-1 1; 3  1];
%A2(:,:,4) = [0 0;  1 -1];
%A2 = A2/4;
V3 = zeros(2,2,3); V3(:,:,2) = V(:,:,1);
Ahat1 = PolyMatConv(U,ParaHerm(V3));
D = zeros(2,2,1); D(:,:,1) = diag([1,-1]);
Ahat2 = PolyMatConv(U,PolyMatConv(D,ParaHerm(V3)));

%------------------------------------------------------------------------------
%  Figure 4
%------------------------------------------------------------------------------
figure(1);
plot((0:1023)/1024,real(fft(s1,1024)),'b-');
hold on;
plot((0:1023)/1024,abs(real(fft(s2,1024))),'r--');
axis([0 1 0 1.5]); grid on;
ylabel('$|\sigma_m(\mathrm{e}^{\mathrm{j}\Omega})|$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/6:1),'XTickLabel',{'$0$','$\pi/3$','$2\pi/3$','$\pi$','$4\pi/3$','$5\pi/3$',...
      '$2\pi$'},...
    'YTick',(0:.5:1.5),'YTickLabel',...
     {'$0$','$.5$','$1$','$1.5$'});
legend({'$m=1$','$m=2$'},'interpreter','latex','fontsize',FS-2,...
        'location','SouthWest');
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
set(gcf,'OuterPosition',[230 250 570 250]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig4.eps

%------------------------------------------------------------------------------
%  Procrustes attempt --- Fourier series approximation
%------------------------------------------------------------------------------
Nfft = 2^8;
IndexN = (-Nfft/2:(Nfft/2-1));
Vprime = zeros(2,2,Nfft);
Vprime(:,1,Nfft/2+1) = [1; -1]/sqrt(2);
% Fourier series of square wave
Vprime(:,2,1:Nfft/2) = [1; 1] * sqrt(2)./(IndexN(1:Nfft/2)*pi).*sin(IndexN(1:Nfft/2)*pi*2/3);
Vprime(:,2,Nfft/2+2:end) = [1; 1] * sqrt(2)./(IndexN(Nfft/2+2:end)*pi).*sin(IndexN(Nfft/2+2:end)*2*pi/3);
Vprime(:,2,Nfft/2+1) = [1; 1]*(1/3)/sqrt(2);
 
%------------------------------------------------------------------------------
%  Figure 5
%------------------------------------------------------------------------------
Qhat3 = PolyMatConv(U,ParaHerm(Vprime));
Qhat3 = real(Qhat3(:,:,117:142));
Anew = zeros(2,2,26);
Anew(:,:,11:14) = A;

figure(5);
t = (-10:15); FS= FS-2;
subplot(221);
stem(t,squeeze(Anew(1,1,:)),'b'); hold on;
plot(t,squeeze(Qhat3(1,1,:)),'r*');
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
ylabel('$a_{1,m}[n], \; \hat{q}^{(N)}_{1,m}[n]$','interpreter','latex','fontsize',FS);
text(-9,0.75,'$m=1$','interpreter','latex','fontsize',FS);
subplot(222);
stem(t,squeeze(Anew(1,2,:)),'b'); hold on;
plot(t,squeeze(Qhat3(1,2,:)),'r*');
text(-9,0.75,'$m=2$','interpreter','latex','fontsize',FS);
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
subplot(223);
stem(t,squeeze(Anew(2,1,:)),'b'); hold on;
plot(t,squeeze(Qhat3(2,1,:)),'r*');
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
text(-9,0.75,'$m=1$','interpreter','latex','fontsize',FS);
ylabel('$a_{2,m}[n], \; \hat{q}^{(N)}_{2,m}[n]$','interpreter','latex','fontsize',FS);
xlabel('time index $n$','interpreter','latex','fontsize',FS);
subplot(224);
stem(t,squeeze(Anew(2,2,:)),'b'); hold on;
plot(t,squeeze(Qhat3(2,2,:)),'r*');
axis([-10.5 10.5 -1.1 1.1]); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-10:5:10),'XTickLabel',{'$-10$','$-5$','$0$','$5$','$10$'});
text(-9,0.75,'$m=2$','interpreter','latex','fontsize',FS);
xlabel('time index $n$','interpreter','latex','fontsize',FS);
set(gcf,'OuterPosition',[230 250 570 350]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig5.eps

error3 = sum(sum(sum(abs(Anew-Qhat3).^2)));

%------------------------------------------------------------------------------
%  Figure 6 --- check paraunitarity via QQ^P
%------------------------------------------------------------------------------
QQ = PolyMatConv(Qhat3,ParaHerm(Qhat3));
n = (-20:20);
QQ = QQ(:,:,6:46);
Ideal = zeros(size(QQ));
Ideal(:,:,21) = eye(2);
AxisDims = [-20 20 -.1 1];
figure(6);
stem(n,squeeze(Ideal(1,1,:)),'bo'); 
hold on; plot(n,squeeze(QQ(1,1,:)),'r*');
axis(AxisDims); grid on;
set(gca,'TickLabelInterpreter','latex','XTick',(-20:5:20),'XTickLabel',{'$-20$','$-10$','$0$','$10$','$20$'});
xlabel('time index $n$','interpreter','latex','fontsize',FS);
ylabel('$r_{Q,1,1}[n]$','interpreter','latex','fontsize',FS);set(gcf,'OuterPosition',[230 250 570 300]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig6.eps

%------------------------------------------------------------------------------
%  error metrics
%------------------------------------------------------------------------------
PowTwos = (3:18);
ErrorNorm = zeros(1,length(PowTwos));
PUError = zeros(1,length(PowTwos));
for i = 1:length(PowTwos),
  Nfft = 2^PowTwos(i);
  IndexN = (-Nfft/2:(Nfft/2-1));
  Vprime = zeros(2,2,Nfft);
  Vprime(:,1,Nfft/2+1) = [1; -1]/sqrt(2);
  % Fourier series of square wave
  Vprime(:,2,1:Nfft/2) = [1; 1] * sqrt(2)./(IndexN(1:Nfft/2)*pi).*sin(IndexN(1:Nfft/2)*pi*2/3);
  Vprime(:,2,Nfft/2+2:end) = [1; 1] * sqrt(2)./(IndexN(Nfft/2+2:end)*pi).*sin(IndexN(Nfft/2+2:end)*2*pi/3);
  Vprime(:,2,Nfft/2+1) = [1; 1]*(1/3)/sqrt(2);

  % solution in the time domain
  Qhat3 = PolyMatConv(U,ParaHerm(Vprime)); 
  % mismatch w.r.t. A
  Error = Qhat3;
  t = (Nfft/2-1:Nfft/2+2);
  Error(:,:,t) = Error(:,:,t) - A;
%  ErrorNorm(i) = sum(sum(sum(abs(Error))));
  for n = 1:size(Error,3),
     %
     %   changed |.|_F to |.|^2_F, SW, 5/1/24
     % 
     ErrorNorm(i) = ErrorNorm(i) + norm(Error(:,:,n),'fro')^2;
  end;   
  % mismatch w.r.t. paraunitarity
  Qf = fft(Qhat3,2*Nfft,3);
  for k = 1:2*Nfft,
     PUError(i) = PUError(i) + sum(sum(abs(squeeze(Qf(:,:,k))*squeeze(Qf(:,:,k))'-eye(2)).^2));
  end;   
  PUError(i) = PUError(i)/(2*Nfft);
end;

%------------------------------------------------------------------------------
%  Figure 7 --- metrics vs order
%------------------------------------------------------------------------------
FS = FS+1;
fig = figure(7); clf;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% for legend only
plot([-1 -1],[1 2],'b-'); hold on; plot([-1 -1],[1 2],'r--'); 
%h = plot([-1 -1],[1 2],'-','linewidth',3);
%set(h(1),'color',[1 1 1]*0.7); 
h = plot([2 18],[1 1]*(17/12-2*sqrt(3)/pi),'-','linewidth',2);
set(h(1),'color',[1 1 1]*0.7); 
plot(PowTwos,ErrorNorm,'r--');
axis([2 18 .2 .325]);
legend({'$e_{\mathrm{PU}}$','$e_{\mathrm{LS}}$',...
    '$\mathrm{lim}_{N\rightarrow \infty}(e_{\mathrm{LS}})$'},'interpreter',...
    'latex','fontsize',FS-2,...
        'location','East');
set(gca,'TickLabelInterpreter','latex','YTick',(0.2:0.025:0.325),'YTickLabel',{'$0.200$','$0.225$','$0.250$','$0.275$','$0.300$','$0.325$'});
%ylabel('$\sum_{n}|\mathbf{A}[n]-\hat{\mathbf{Q}}_[n]|^2_{\mathrm{F}}$',...
%       'interpreter','latex','fontsize',FS);
ylabel('least squares mismatch, $e_{\mathrm{LS}}$',...
       'interpreter','latex','fontsize',FS);
grid on;
yyaxis right
plot(PowTwos,10*log10(PUError),'b-');
axis([2 18 -55 -5]);
set(gca,'TickLabelInterpreter','latex','YTick',(-55:10:-5),'YTickLabel',{'$-55$','$-45$','$-35$','$-25$','$-15$','$-5$'});
%ylabel('$10\mathrm{log}_{10}\{\sum_n \|\mathbf{Q}^\prime[n]\ast\mathbf{Q}^{\prime\mathrm{H}}[-n]-\mathbf{I}\delta[n]%%\|^2_{\mathrm{F}}\}$','interpreter','latex','fontsize',FS);
ylabel('paraunitarity error, $10\mathrm{log}_{10}\{ e_{\mathrm{PU}} \}$','interpreter','latex','fontsize',FS);
xlabel('$\mathrm{log}_2 N$','interpreter','latex','fontsize',FS);
set(gcf,'OuterPosition',[230 250 570 350]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig7.eps

