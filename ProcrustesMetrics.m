function Results = ProcrustesMetrics(A,Sigma,U,V,Qhat);
% ProcrustesMetrics(A,Sigma,U,V,Qhat);
%
% This function evaluates a few metrics for the Procrustes solution Qhat to a 
% matrix of analytic functions A, whose analytic SVD U(z) Sigma(z) V^P(z) is
% given.
%
% Input parameters:
%    A           matrix of impulse responses / polynomials
%    Sigma       analytic singular values
%    U           analytic left-singular vectors
%    V           analytic right-singular vectors
%    Qhat        approximated Procrustes solution
%
% Output parameter:
%    Results     row vectors with errors in paraunitarity, diagonalisation, 
%                  positivity, and the least squares mismatch between A and 
%                  Qhat.
 
% S. Weiss, UoS, 14/6/2024 


%-----------------------------------------------------------------------------
%  paraunitarity and diagonalisation
%-----------------------------------------------------------------------------
% metric --- paraunitarity
QQ = PolyMatConv(Qhat,ParaHerm(Qhat));
Lq = size(QQ,3);
QQ(:,:,(Lq+1)/2) = QQ(:,:,(Lq+1)/2) - eye(3);
paraunitarity = PolyMatNorm(QQ);
disp(sprintf('paraunitarity error        %0.5g',paraunitarity));

% metric --- diagonality
Pi = PolyMatConv(ParaHerm(U),PolyMatConv(Qhat,V));            %  extract switching functions
S = PolyMatConv(Sigma,Pi);
diagonality = PolyMatNorm(S,'OffDiag')/PolyMatNorm(S,'OnDiag');
disp(sprintf('error in diagonalisation   %0.5g',diagonality));

% metric --- positivity
% S is mostly non-negative and real on the unit circle; hence in the time 
% domain is will satisfy approximately the properties of an autocorrelation
% sequence. Hence we look for its maximum to determine zero lag.
M = size(A,1);
LS = size(S,3);
dummy = zeros(LS,1);
NDFT = 2^(ceil(log2(LS))+2); 
SS = zeros(NDFT,M);
for m = 1:M,
   dummy2 = squeeze(S(m,m,:));   
   SS(1:LS,m) = dummy2;
   dummy = dummy + abs(dummy2);
end;
[~,MaxIndex] = max(dummy);
Advance = MaxIndex-1;
disp(sprintf('required advance           %d',Advance));
SS = circshift(SS,-Advance,1);
SSf = fft(SS,NDFT,1);
int_S = sum(sum(real(SSf).^2));
int_Sigma = 0;
for m = 1:M,
    int_Sigma = int_Sigma + sum(abs(fft(squeeze(Sigma(m,m,:)),NDFT).^2));
end;
positivity = (int_Sigma - int_S)/int_Sigma;    
disp(sprintf('positivity error           %0.5g',positivity));

if nargout==3,
   Results = [paraunitarity diagonality positivity];
else
   [~,dummy,~,~] = PolyMatAlign(A,Qhat);
   Results = [paraunitarity diagonality positivity dummy];
end;   
