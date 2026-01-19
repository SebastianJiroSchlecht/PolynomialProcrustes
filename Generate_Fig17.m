% Generate Figure 17 for the paper
%    S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary 
%    Approximation of Matrices of Analytic Functions," submitted to IEEE 
%    Trans. Signal Process., subitted Mar. 2025
%
% show singular values for a matrix with M=32

%------------------------------------------------------------------------------
%  parameters
%------------------------------------------------------------------------------
Nfft = 1024;
TrimThreshold = 1e-10;
SeedVal = 1;
DispMode='on';
M = 32;
L = 1;

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
dummy = zeros(M,Nfft); Ls = size(s,2);
dummy(:,1:(Ls+1)/2) = s(:,(Ls+1)/2:Ls);
dummy(:,end-(Ls-1)/2+1:end) = s(:,1:(Ls+1)/2-1);  
sf =  real(fft(dummy,Nfft,2));
if strcmp(DispMode,'on')==1,
   figure(1); clf;
   plot((0:Nfft-1)/Nfft,sf');
   xlabel('norm. freq.'); ylabel('singular values'); 
end;   
S = zeros(M,M,size(s,2));
for m = 1:M;
  S(m,m,:) = s(m,:);
end;  

%------------------------------------------------------------------------------
%  Figure 17 --- plot singular values
%------------------------------------------------------------------------------
figure(17);
st = zeros(M,1024);
st(:,1:2) = s(:,2:3); st(:,1024)=s(:,1);
sf = real(fft(st,1024,2));
plot((0:1023)/1024,sf);
axis([0 1 -3.5 2.25]);
set(gca,'TickLabelInterpreter','latex','XTick',(0:1/8:1),'XTickLabel',...
   {'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$',...
   '$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'});
xlabel('normalised angular frequency $\Omega$','interpreter','latex'); 
ylabel('$\sigma_m(\mathrm{e}^{\mathrm{j}\Omega})$','interpreter','latex');
grid;
set(gcf,'OuterPosition',[230 250 570 300]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc TSP_PPP_Fig17.eps


