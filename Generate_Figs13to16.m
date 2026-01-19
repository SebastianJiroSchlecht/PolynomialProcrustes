% Generate Figures 13, 14, 15, and 16 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., submitted Mar. 2025
%
% simulation for one specific matrix
 
clear all; close all;

M = 3; N = 10; FS = 12;

%------------------------------------------------------------------------------
%  matrix A
%------------------------------------------------------------------------------
% singular values
Nfft = 64;
s = [0  1 0  1 0;
     0 .5 0 .5 0; 
     -1i*.5 .5 0 .5 1i*.5];
st = zeros(3,Nfft); st(:,1:3) = s(:,3:5); st(:,Nfft-1:Nfft) = s(:,1:2);     
sf = fft(st,Nfft,2);
figure(1); plot((0:Nfft-1)/Nfft,real(sf)','*-'); 

U = PUPolyMatRand(3,10,0,'complex');
V = PUPolyMatRand(3,10,1,'complex');
S = zeros(3,3,5);
for m = 1:3, S(m,m,:) = s(m,:); end;
A2 = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

B = zeros(M,M,1); B(:,:,1) = eye(M);
A = zeros(M,M,size(A2,3)+N+1);
A(:,:,N+2:end) = A2;
if exist('Qhat1024.mat','file'),
   load Qhat1024.mat
else    
   Qhat = PUProcrustes(A,B,2^12,0,15);
   save Qhat1024.mat Qhat             % Procrustes 
end;

%-----------------------------------------------------------------------------
%  Figure 13: singular values
%-----------------------------------------------------------------------------
figure(13); clf; 
Nfft2 = 1024; Nfft3 = 16; 
ss = zeros(3,Nfft2); 
ss(:,1:3) = s(:,3:5); ss(:,Nfft2-1:Nfft2) = s(:,1:2); 
Sf = real(fft(ss,Nfft2,2));
Sf = [Sf Sf(:,1)];
% plot curves
plot((0:Nfft2)/Nfft2,Sf(1,:),'b-');
hold on;
plot((0:Nfft2)/Nfft2,Sf(2,:),'r--');
h = plot((0:Nfft2)/Nfft2,Sf(3,:),'-.');
set(h(1),'color',[0 0.5 0]);
plot((0:Nfft2/Nfft3:Nfft2)/Nfft2,Sf(1,1:Nfft2/Nfft3:end),'b*');
plot((0:Nfft2/Nfft3:Nfft2)/Nfft2,Sf(2,1:Nfft2/Nfft3:end),'r*');
h = plot((0:Nfft2/Nfft3:Nfft2)/Nfft2,Sf(3,1:Nfft2/Nfft3:end),'*');
set(h(1),'color',[0 0.5 0]);
axis([0 1 -2 2]); grid on;
plot([0 1 2 3 4]/4,[1 0 -1 0 1],'ko','MarkerFaceColor','k');
plot([0.08333449 .5-0.08333449],[1 -1]*1.732036,'ko','linewidth',2','MarkerFaceColor','w');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-2:1:2),'YTickLabel',...
     {'$-2$','$-1$','$0$','$1$','$2$'});
xlabel('normalised anguar frequency $\Omega$','interpreter','latex','fontsize',FS);
ylabel('$\sigma_{m}(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',FS);
legend({'$m=1$','$m=2$','$m=3$'},...
       'interpreter','latex','fontsize',FS-2,'location','SouthWest');	
set(gcf,'OuterPosition',[230 250 570 280]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig13.eps

Qhat = Qhat(:,:,513:end);
%-----------------------------------------------------------------------------
%  Figure 14: matrix A(z)
%-----------------------------------------------------------------------------
figure(14); clf;
for m = 1:3,
   for n = 1:3,
      subplot(3,3,3*(m-1)+n); 
      stem((0:24),squeeze(real(A2(m,n,:))),'b'); hold on;
      plot((0:24),squeeze(real(Qhat(m,n,12:36))),'r*'); 
      axis([0 24 -.4 .4]); grid on;
      if n == 1,
        ylabel(sprintf('$\\ell=%d$',m),'interpreter','latex'); 
      end;   
      if m==1,
         title(sprintf('$m=%d $',n),'interpreter','latex');
      end;
      if m==3,
         xlabel('index $n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex',...
      'XTick',(0:8:24),'XTickLabel',{'$0$','$8$','$16$','$24$'},...
      'YTick',(-.4:.2:.4),'YTickLabel',...
     {'$-.4$','$-.2$','$0$','$.2$','$4$'});
   end;
end;
set(gcf,'OuterPosition',[230 250 570 400]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig14.eps

%-----------------------------------------------------------------------------
%  Figure 16: \Pi(z) \Sigma(z)
%-----------------------------------------------------------------------------
P = PolyMatConv(ParaHerm(U),PolyMatConv(Qhat,V));
Ls = size(S,3);
Pi = zeros(1024,3);
% need to account for non-causal allpass
for i = 1:3,
   dummy = zeros(1024,1); 
   dummy(1:size(P,3)) = squeeze(P(i,i,:));
   Pi(:,i) = circshift(dummy,-33);  
%   Pi(1024-Shift:end,i) = squeeze(P(i,i,1:Shift)); 
end;
Pif = fft(Pi,1024,1);
Pif = [Pif; Pif(1,:)];
SigmaPif = real(Pif.*Sf.');

% metric --- diagonalisation
SSS = PolyMatConv(S,P);
diagon = PolyMatNorm(SSS,'OffDiag')/PolyMatNorm(SSS,'OnDiag');
% metric --- positivity
positivity = norm(SigmaPif-abs(Sf.'),'fro').^2/norm(abs(Sf),'fro').^2;
% metrix --- paraunitarity
QQ = PolyMatConv(Qhat,ParaHerm(Qhat));
Lq = size(QQ,3);
QQ(:,:,(Lq+1)/2) = QQ(:,:,(Lq+1)/2) - eye(3);
paraunitarity = PolyMatNorm(QQ);

figure(16); clf;
% plot for legend
plot([-1 -2],[0 0],'b-'); hold on; 
plot([-1 -2],[0 0],'r--');  
h = plot([-1 -2],[0 0],'-.');  set(h(1),'color',[0 0.5 0]);
% plot curves
for i = 1:3,
   h = plot((0:Nfft2)/Nfft2,abs(Sf(i,:)),'-','linewidth',3); 
   set(h(1),'color',[1 1 1]*0.7); 
end;   
plot((0:Nfft2)/Nfft2,SigmaPif(:,1),'b-');
hold on;
plot((0:Nfft2)/Nfft2,SigmaPif(:,2),'r--');
h = plot((0:Nfft2)/Nfft2,SigmaPif(:,3),'-.');
set(h(1),'color',[0 0.5 0]);
axis([0 1 0 2]); grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(0:1:2),'YTickLabel',...
     {'$0$','$1$','$2$'});
xlabel('normalised anguar frequency $\Omega$','interpreter','latex','fontsize',FS);
ylabel('$|\sigma_{m}(\mathrm{e}^{\mathrm{j}\Omega})|$, $\Re\{s_{m,m}(\mathrm{e}^{\mathrm{j}\Omega})\}$',...
	'interpreter','latex','fontsize',FS);
legend({'$m=1$','$m=2$','$m=3$'},...
       'interpreter','latex','fontsize',FS-2,'location','NorthEast');	
Box1 = [.245 .255 0 0.05];
plot(Box1([1 1 2 2 1]),Box1([3 4 4 3 3]),'k-','linewidth',1);
axes('position',[.29 .64 0.07 .25]);
for i = 1:3,
   h = plot((0:Nfft2)/Nfft2,abs(Sf(i,:)),'-','linewidth',2); 
   set(h(1),'color',[1 1 1]*0.75); hold on;
end;  
plot((0:Nfft2)/Nfft2,SigmaPif(:,1),'b-');
hold on;
plot((0:Nfft2)/Nfft2,SigmaPif(:,2),'r--');
h = plot((0:Nfft2)/Nfft2,SigmaPif(:,3),'-.');
set(h(1),'color',[0 0.5 0]); 
axis([0.245 0.255 0 0.05]); grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',[0.245 0.255],'XTickLabel',{'$\frac{49\pi}{100}$','$\frac{51\pi}{100}$'},...
    'YTick',[0 0.05],'YTickLabel',...
     {'$0.00$','$0.05$'});
set(gcf,'OuterPosition',[230 250 570 280]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig16.eps

%-----------------------------------------------------------------------------
%  Figure 15: polynomial Procrustes solution
%-----------------------------------------------------------------------------
figure(15); clf;
for m = 1:3,
   for n = 1:3,
      subplot(3,3,3*(m-1)+n); 
      stem((0:24),squeeze(imag(A2(m,n,:))),'b'); hold on;
      plot((0:24),squeeze(imag(Qhat(m,n,12:36))),'r*'); hold on;
      axis([0 24 -.4 .4]); grid on;
      if n == 1,
%         ylabel(sprintf('$\\hat{q}_{\\ast,%d,m}[n]$',m),'interpreter','latex'); 
         ylabel(sprintf('$\\ell=%d$',m),'interpreter','latex'); 
      end;   
      if m==1,
         title(sprintf('$m=%d $',n),'interpreter','latex');
      end;
      if m==3,
         xlabel('index $n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex',...
      'XTick',(0:8:24),'XTickLabel',{'$0$','$8$','$16$','$24$'},...
      'YTick',(-.4:.2:.4),'YTickLabel',...
     {'$-.4$','$-.2$','$0$','$.2$','$4$'});
   end;
end;
set(gcf,'OuterPosition',[230 250 570 400]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig15.eps  
      
%-----------------------------------------------------------------------------
%  some numerical evaluations
%-----------------------------------------------------------------------------
Res = ProcrustesMetrics(A,S,U,V,Qhat);   
disp(sprintf('error in paraunitarity:     %2.12g',Res(1)));
disp(sprintf('error in diagonalisation:   %2.12g',Res(2)));
disp(sprintf('error in positivity:        %2.12g',Res(3)));
disp(sprintf('least squares error A-Q:    %2.12g',Res(4)));
[~,dummy,~,~] = PolyMatAlign(A,PolyMatConv(U,ParaHerm(V)));
disp(sprintf('least squares error A-UV^P: %2.12g',dummy));

