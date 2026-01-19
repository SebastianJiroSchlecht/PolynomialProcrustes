% PUProcrustesTest.m
%
% testbench for function PUProcrustes()

% S. Weiss, 29/3/2025

TEST = 1;

switch TEST,
   case 0,
      disp('test case 1: random 2x2 matrix');
      randn('seed',0);
      A = randn(2,2,2);
      B = zeros(2,2,1); B(:,:,1) = eye(2);
      Q = PUProcrustes(A,B,1024,0,20);      
   case 1, 
      disp('test case 2: oversampled matrix, which remedies a single singular values with a single zero crossing');    
      A = zeros(1,1,3);
      A(1,1,:) = [1 0 1]; 
      B = zeros(1,1,1); B(1,1,1) = 1;
      Q = PUProcrustes(A,B,1024,0,20);  
   case 2, 
      disp('test case 2: diagonal 2x2 matrix, singular values with zero crossings');    
      A = zeros(2,2,3);
      A(1,1,:) = [1 0 1]; A(2,2,:) = [0 1 0];
      B = zeros(2,2,1); B(:,:,1) = eye(2);
      Q = PUProcrustes(A,B,1024,0,20);  
   otherwise,
      disp('test case not implemented');
end;    
