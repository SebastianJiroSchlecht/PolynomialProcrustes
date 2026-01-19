function Results = ProcrustesRandomMatrixSingleTest(SeedVal);
% ProcrustesRandomMatrixSingleTest.m
% 
% Perform one simulation run within an ensemble of tests using randomised
% matrices with ground truth Procrustes solution. The randomisation is 
% initialised by the seed value SeedVal. The function returns a number of
% metrics for the simulations in [1] in the variable Results:
%    Results(1)       seed value
%    Results(2)       min. LS mismatch
%    Results(3)       zeta from iteration
%    Results(4)       DFT length
%    Results(5)       paraunitarity error
%    Results(6)       LS mismatch
% 
% Input parameter:
%       SeedVal            seed value for random number generator
%
% Output parameter:
%       Results              vector containing various metrics
%
% [1] S. Weiss, S.J. Schlecht, M. Moonen: "Best Least Squares Paraunitary 
%     Approximation of Matrices of Analytic Functions", submitted to IEEE
%     Trans. Signal Process., March 2025.

%-------------------------------------------------------------------
%   generate matrices
%-------------------------------------------------------------------
M = 6;                          %   spatial dimensions
N = 64;                        %   order of allpass design
[A2,Q,U,S,V] = ProcrustesRandomMatrix(M,6,N,SeedVal,'on');

%-------------------------------------------------------------------
%   calculate minimum least squares error 
%-------------------------------------------------------------------
LS3 = size(S,3);
SS = zeros(M,2^14);
for m = 1:M,
   SS(m,1:LS3) = permute(S(m,m,:),[1 3 2]);
end;
SSf = real(fft(circshift(SS,[0 -(LS3-1)/2]),2^14,2));
% numerical integration
MinLSError = (norm(abs(SSf)-ones(M,2^14),'fro')^2 )/(2^14);
disp(sprintf('min. least squares mismatch     %0.5g',MinLSError));
Results = zeros(1,6);
Results(1) = SeedVal;
Results(2) = MinLSError;   

%-------------------------------------------------------------------
%   determine polynomial Procrustes solution 
%-------------------------------------------------------------------
B = zeros(M,M,1); B(:,:,1) = eye(M);
A = zeros(M,M,size(A2,3)+2*N+1);
A(:,:,2*N+2:end) = A2;    % delay to enable a  causal solution      
[Qhat,~,~,~,Results(3),Results(4)] = PUProcrustes(A,B,2^14,0,N);      % Procrustes 
PMetrics = ProcrustesMetrics(A,S,U,V,Qhat);
disp(sprintf('paraunitarity error             %0.5g',PMetric(1)));
disp(sprintf('least squares mismatch          %0.5g',PMetric(4)));
Results(5) = PMetric(1); Results(6) = PMetric(4);
