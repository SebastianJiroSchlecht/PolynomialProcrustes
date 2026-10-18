function test_reproducibility
% Small numerical regressions; does not regenerate the Monte Carlo ensembles.
addpath(fullfile(pwd,'PEVDToolbox','General'));
addpath(fullfile(pwd,'PEVDToolbox','Decompositions'));
set(groot,'defaultFigureVisible','off');
cleanup = onCleanup(@() close('all'));

% Metrics must work for dimensions used outside the three-channel example.
for M = [2 3 6]
    I = eye(M);
    r = ProcrustesMetrics(I,I,I,I,I);
    assert(all(abs(r)<1e-12),'Identity metrics must vanish for every dimension.');
end

% Defaults and rectangular A/B must agree with explicit normal-equation input.
B = diag([2 1]);
Q = PUProcrustes(B);
QQ = PolyMatConv(Q,ParaHerm(Q));
QQ(:,:,(size(QQ,3)+1)/2) = QQ(:,:,(size(QQ,3)+1)/2)-eye(2);
assert(PolyMatNorm(QQ)<1e-10);
A = [eye(2); zeros(1,2)];
B = [diag([2 1]); zeros(1,2)];
Qrect = PUProcrustes(B,A,64,0,32);
Qsquare = PUProcrustes(A'*B,eye(2),64,0,32);
[~,err] = PolyMatAlign(Qrect,Qsquare);
assert(err<1e-12,'Rectangular normal equations must match square input.');

% Missing bins include repeated singular values, as well as spectral zeros.
At = zeros(2,2,3);
At(1,1,:) = [.1 2 .1]; At(2,2,:) = [0 2 0];
[Q,~,~,flag] = PUProcrustes(At,eye(2),64,1,8);
assert(flag==1 && all(isfinite(Q(:))));

% This used to fail on a misspelled variable after computing the metrics.
Qloose = PUProcrustes(diag([2 1]),eye(2),64,0,8,1e-4);
assert(all(isfinite(Qloose(:))));
r = ProcrustesRandomMatrixSingleTest(0);
assert(numel(r)==6 && all(isfinite(r)));
assert(r(3)>=0 && r(3)<1 && r(4)>=16 && mod(log2(r(4)),1)==0,...
    'Columns 3 and 4 must contain convergence error and FFT length.');
fprintf('Numerical regressions passed.\n');
end
