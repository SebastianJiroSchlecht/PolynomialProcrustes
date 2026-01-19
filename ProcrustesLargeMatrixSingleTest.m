function Results = ProcrustesLargeMatrixSingleTest(SeedVal,M);
% ProcrustesLargeMatrixSingleTest.m
% 
% Perform one simulation run within an ensemble of tests using randomised
% matrices with ground truth Procrustes solution. The randomisation is 
% initialised by the seed value SeedVal. The function returns a number of
% metrics for the simulations in [1] in the variable Results:
%    Results(1)       seed value
%    Results(2)       M
%    Results(3)       DFT length
%    Results(4)       length of paraunitary matrix
%    Results(5)       paraunitarity error
%    Results(6)       LS mismatch
%    Results(7)       computation time 
%
% Input parameter:
%       SeedVal            seed value for random number generator
%       M                  spatial dimension of matrix
%
% Output parameter:
%       Results              vector containing various metrics
%
% [1] S. Weiss, S.J. Schlecht, M. Moonen: "Best Least Squares Paraunitary 
%     Approximation of Matrices of Analytic Functions", submitted to IEEE
%     Trans. Signal Process., March 2025.

%-------------------------------------------------------------------
%   parameters
%-------------------------------------------------------------------
N = 64;                        %   order of allpass design
Nfft = 1024*2^5;
L = 1;
DispMode = 'off';

%------------------------------------------------------------------------------
%  generate analytic SVD factors
%------------------------------------------------------------------------------
U = PUPolyMatRand(M,L,SeedVal,'complex');
V = PUPolyMatRand(M,L,SeedVal+1,'complex');
randn('seed',SeedVal);
% order of singular values must be even
s = randn(M,ceil(L/2)*2+1) + 1i*randn(M,ceil(L/2)*2+1);
s = (s + fliplr(conj(s)))/2;
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
%  Procrustes solution
%------------------------------------------------------------------------------
B = zeros(M,M,1); 
B(:,:,1) = eye(M);
tstart = tic;
[Q,D,Sfcorr2,OSflag,Nfft,e] = PUProcrustes(A,B,Nfft,0,16);
CompTime = toc(tstart);
%------------------------------------------------------------------------------
%  truncate outer zeros of paraunitary matrix
%------------------------------------------------------------------------------
MaxElement = zeros(size(Q,3),1);
% element size
for i = 1:size(Q,3),
   MaxElement(i) = max(max(abs(Q(:,:,i))));
end;   
% find any leading small components
StartIndex = 1; EndIndex = size(Q,3);
while MaxElement(StartIndex)<1e-8, StartIndex = StartIndex+1; end;
while MaxElement(EndIndex)<1e-8, EndIndex = EndIndex-1; end;
Q = Q(:,:,StartIndex:EndIndex);  

%------------------------------------------------------------------------------
%  metrics
%------------------------------------------------------------------------------
% paraunitarity
QQ = PolyMatConv(Q,ParaHerm(Q));
Lq = size(QQ,3);
QQ(:,:,(Lq+1)/2) = QQ(:,:,(Lq+1)/2) - eye(M);
paraunitarity = PolyMatNorm(QQ);
disp(sprintf('paraunitarity error        %0.5g',paraunitarity));

% mismatch
[Shift,error,A2,B2] = PolyMatAlign(A,Q);

%------------------------------------------------------------------------------
%  assign outputs
%------------------------------------------------------------------------------
Results(1) = SeedVal;
Results(2) = M;
Results(3) = Nfft;
Results(4) = size(Q,3);
Results(5) = paraunitarity;
Results(6) = error;
Results(7) = CompTime;

