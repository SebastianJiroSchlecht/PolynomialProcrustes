function [A,Q,U,S,V] = ProcrustesRandomMatrix(M,L,N,SeedVal,DispMode);
% ProcrustesRandomTestMatrix(M,L,N,SeedVal,DispMode);
% 
% [A,Q,S,U,V] = ProcrustesRandomMatrix(M,L,SeedVal) generates a random MxM
% matrix A of transfer functions, with ground truth analytic SVD factors U, S,
% and V all of order L. These are randomised functions initialised by the seed 
% value SeedVal.
%
% Q is the polynomial Procrustes solution with switching functions of order N,
% based on the knowledge of the ground truth matrices U,S, and V. This solution
% is truncated to a coefficient threshold that is internal to the function.
%
% Input parameters
%    M          spatial dimension of matrices A, Q, U, S, V
%    L          orders of U, S, and V
%    N          order of allpass filters
%    SeedVal    initialisation for random number generator
%    DispMode   display mode {'on','off'} 
%                           (optional; default is 'off')
%
% Output parameters
%    A          MxMx(3L+1) matrix of transfer functions
%    Q          MxMxLq polynomial Procrustes solution to A
%    U,S,V      MxMx(L+1) analytic singular values and vectors of order L

% S. Weiss, UoS, 3/3/24

%------------------------------------------------------------------------------
%  parameters and check
%------------------------------------------------------------------------------
if nargin<5,
   DispMode='off';
end;
Nfft = 1024;
TrimThreshold = 1e-10;

%------------------------------------------------------------------------------
%  generate analytic SVD factors
%------------------------------------------------------------------------------
U = PUPolyMatRand(M,L,SeedVal,'complex');
V = PUPolyMatRand(M,L,SeedVal+1,'complex');
randn('seed',SeedVal);
% order of singular values must be even
s = randn(M,floor(L/2)*2+1) + 1i*randn(M,floor(L/2)*2+1);
s = s + fliplr(conj(s));
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
%  matrix of transfer functions
%------------------------------------------------------------------------------
A = PolyMatConv(U,PolyMatConv(S,ParaHerm(V)));

%------------------------------------------------------------------------------
%  optimum Procrustes solution 
%------------------------------------------------------------------------------
% determine switching functions
sfs = sign(sf);
At = zeros(M,Nfft);
for m = 1:M,
   i = 0;
   % find all zero crossings in m-th singular value
  Omega = [];
   for k = 2:Nfft,
       if sfs(m,k)~=sfs(m,k-1),
          i = i+1;
          Omega(i) = ZeroCrossing(s(m,:),2*pi*(k-2)/Nfft,2*pi*(k-1)/Nfft);
       end;
       % if there is a zero crossing between the last and the first bin, we have missed it!
       %     to address this:
      if sfs(m,1)~=sfs( m,Nfft),
          i = i+1;
        Omega(i) =  ZeroCrossing(s(m,:),2*pi*(Nfft-1)/Nfft,2*pi);
      end; 
   end;
   % display if asked
   if strcmp(DispMode,'on')==1,
      figure(1); hold on;
      plot(Omega/2/pi,zeros(size(Omega)),'ko');
   end;
   if length(Omega)>0,
      [b,a] = SingularValueAllpassSwitchCompact(Omega,N);  
      Delays(m) = N-length(Omega)/2;            
%      [b,a] = SingularValueAllpassSwitch(Omega,N);  
 %     Delays(m) = (N-1)*length(Omega)/2;    
       At(m,:) = impz(b,a,Nfft).';         % non-causal allpass filter coefficients
   else                                   % singular values has NO zero crossings
      At(m,1) = 1;
      Delays(m) = 0;
   end;
   Af(m,:) = fft(circshift(At(m,:),-Delays(m),2),Nfft,2);  
end;
% assemble Procrustes solution in the DFT domain
Uf = fft(U,Nfft,3); Vf = fft(V,Nfft,3);
Qf = zeros(M,M,Nfft);
for k = 1:Nfft,
   Qf(:,:,k) = Uf(:,:,k)*diag(Af(:,k))*Vf(:,:,k)';
end;   
% convert to time domain and adjust
Q = ifft(Qf,Nfft,3);
Q = circshift(Q,max(Delays)+N,3); 
Q = PUPolyMatTrim(Q,TrimThreshold);


%------------------------------------------------------------------------------
%  function ZeroCrossing() 
%------------------------------------------------------------------------------
function Omega = ZeroCrossing(s,Omega1,Omega2);
% Refine the localisation of a zero crossing between the frequencies Omega1
% and Omega2, for a function with a real-values spectrum given by its time-
% domain coefficients s. 
Ls = length(s);
nn = (-(Ls-1)/2:(Ls-1)/2)';
S1 = real(s*exp(-1i*Omega1*nn));
S2 = real(s*exp(-1i*Omega2*nn));
% bisection method for refinement of zero location
if S2<S1,             % function is falling; convert to rising
   s = -s; 
end;   
for i = 1:10,
   Omega = (Omega1+Omega2)/2;
   S12 = real(s*exp(-1i*Omega*nn));
   if S12<0;
      Omega1 = Omega;
   else 
      Omega2 = Omega;   
   end;
end     
