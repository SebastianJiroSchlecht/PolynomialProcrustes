% ProcrustesLargeMatrixEnsembleTest.m
%
% Matlab script file to generate the ensemble simulation for [1].
%
% [1] S. Weiss, S.J. Schlecht, M. Moonen: "Best Least Squares Paraunitary
%     Approximation of Matrices of Analytic Functions," submitted to IEEE
%     Trans. Signal Process., Mar. 2025.

MVals = [2, 4, 8, 16, 32, 64];
%MVals = [64];
SeedVals = (11:15);
for m = 1:length(MVals),
  disp(sprintf('Matrix dimension: %d',MVals(m)));
  for i = 1:length(SeedVals),
    disp(sprintf('Seedvalue: %d',SeedVals(i)));
    Results = ProcrustesLargeMatrixSingleTest(SeedVals(i),MVals(m));
    
    %--------------------------------------------------------
    %  write results to a file
    %--------------------------------------------------------
    if exist('LargeMatrixEnsembleResults.txt','file') ~= 2,
      disp('new results file created');
      dlmwrite('LargeMatrixEnsembleResults.txt',[Results(1) Results(2) Results(3) Results(4) Results(5) Results(6) Results(7)]);
    else  
      disp('results appended');
      dlmwrite('LargeMatrixEnsembleResults.txt',[Results(1) Results(2) Results(3) Results(4) Results(5) Results(6) Results(7)],'-append');
    end;
  end;  
end;
