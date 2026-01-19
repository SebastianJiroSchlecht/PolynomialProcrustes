% ProcrustesRandomMatrixEnsembleTest.m
%
% Matlab script file to generate the ensemble simulation for [1].
%
% [1] S. Weiss, S.J. Schlecht, M. Moonen: "Best Least Squares Paraunitary
%     Approximation of Matrices of Analytic Functions," submitted to IEEE
%     Trans. Signal Process., Mar. 2025.

SeedVals = (1:1000);
for i = 1:length(SeedVals),
   disp(sprintf('Seedvalue: %d',SeedVals(i)));
   Results = ProcrustesRandomMatrixSingleTest(SeedVals(i));
    
   %--------------------------------------------------------
   %  write results to a file
   %--------------------------------------------------------
   if exist('EnsembleResults.txt','file') ~= 2,
      disp('new results file created');
      dlmwrite('EnsembleResults.txt',[Results(1) Results(2) Results(3) Results(4) Results(5) Results(6)]);
   else  
      disp('results appended');
      dlmwrite('EnsembleResults.txt',[Results(1) Results(2) Results(3) Results(4) Results(5) Results(6)],'-append');
   end;
end;
