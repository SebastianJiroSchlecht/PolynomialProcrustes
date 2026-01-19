% Generate Figures 19 and 20 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., submitted Mar. 2025 
%
% This figure compare a few elements of a large matrix A(z) with its Procrustes 
% solution. 

clear all; close all;
%------------------------------------------------------------------------------
%  parameters
%------------------------------------------------------------------------------
Nfft = 1024*4;
SeedVal = 1;
M = 96;
L = 1;
% rows and columns to plot
Rows=[16 32 55]; Columns = [93 94 95 96];

if (exist('LargeMatrixM96Example.mat')~=2)
  %------------------------------------------------------------------------------
  %  generate analytic SVD factors
  %------------------------------------------------------------------------------
  U = PUPolyMatRand(M,L,SeedVal,'complex');
  V = PUPolyMatRand(M,L,SeedVal+1,'complex');
  randn('seed',SeedVal);
  % order of singular values must be even
  s = randn(M,ceil(L/2)*2+1) + 1i*randn(M,ceil(L/2)*2+1);
  s = (s + fliplr(conj(s)))/4;
  s = diag(sign(sum(s,2)))*s/L;      % ensure that for Omega=0 all SVs are nonnegative
  S = zeros(M,M,size(s,2));
  for m = 1:M;
     S(m,m,:) = s(m,:);
  end;  

  %------------------------------------------------------------------------------
  %  matrix of transfer functions
  %------------------------------------------------------------------------------
  A = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

  %------------------------------------------------------------------------------
  %  Procrustes solution
  %------------------------------------------------------------------------------
  B = zeros(M,M,1); 
  B(:,:,1) = eye(M);
  [Q,D,Sfcorr2,OSflag,Nfft,e] = PUProcrustes(A,B,Nfft,0,64);

  %------------------------------------------------------------------------------
  %  align against original solution
  %------------------------------------------------------------------------------
  [Shift,error,A2,Q2] = PolyMatAlign(A,Q);


  %------------------------------------------------------------------------------
  %   extract some elements of matrices
  %------------------------------------------------------------------------------
  t = (-5:10);
  AA = zeros(M,M,16);
  AA(:,:,6:10) = A(:,:,1:5)/2;
  Shift=2048;
  QQ = Q(:,:,Shift-4:Shift+11);
  save LargeMatrixM96Example.mat AA QQ
else
  load LargeMatrixM96Example.mat
end;    

%------------------------------------------------------------------------------
%  Figure 19 --- compare real parts
%------------------------------------------------------------------------------
figure(19);
AA = AA(Rows,Columns,:); QQ = QQ(Rows,Columns,:);
t = (-5:10);
for m = 1:length(Rows),
   for n = 1:length(Columns),
      subplot(length(Rows),length(Columns),length(Columns)*(m-1)+n); 
      stem(t,squeeze(real(AA(m,n,:))),'b'); hold on;
      plot(t,squeeze(real(QQ(m,n,:))),'r*'); 
      axis([-5 10 -.1 .1]); grid on;
      if n == 1,
        ylabel(sprintf('$\\ell=%d$',Rows(m)),'interpreter','latex'); 
      end;   
      if m==1,
         title(sprintf('$m=%d $',Columns(n)),'interpreter','latex');
      end;
      if m==length(Rows),
         xlabel('index $n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex',...
      'XTick',(-5:5:10),'XTickLabel',{'$-5$','$0$','$5$','$10$'},...
      'YTick',(-.1:.05:.1),'YTickLabel',...
     {'$-.1$','$-.05$','$0$','$0.05$','$0.1$'});
   end;
end;
set(gcf,'OuterPosition',[230 250 570 400]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig19.eps

%------------------------------------------------------------------------------
%  Figure 20 --- compare imaginary parts
%------------------------------------------------------------------------------
figure(20);
for m = 1:length(Rows),
   for n = 1:length(Columns),
      subplot(length(Rows),length(Columns),length(Columns)*(m-1)+n); 
      stem(t,squeeze(imag(AA(m,n,:))),'b'); hold on;
      plot(t,squeeze(imag(QQ(m,n,:))),'r*'); 
      axis([-5 10 -.1 .1]); grid on;
      if n == 1,
        ylabel(sprintf('$\\ell=%d$',Rows(m)),'interpreter','latex'); 
      end;   
      if m==1,
         title(sprintf('$m=%d $',Columns(n)),'interpreter','latex');
      end;
      if m==length(Rows),
         xlabel('index $n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex',...
      'XTick',(-5:5:10),'XTickLabel',{'$-5$','$0$','$5$','$10$'},...
      'YTick',(-.1:.1:.1),'YTickLabel',...
     {'$-.1$','$-.05$','$0$','$0.05$','$0.1$'});
   end;
end;
set(gcf,'OuterPosition',[230 250 570 400]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig20.eps


