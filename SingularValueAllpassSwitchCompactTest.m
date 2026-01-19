% SingularValueAllpassSwitchCompactTest.m
%
clear all; close all;

Nfft = 1024;
sigma = [-1i/2 1/4 1/2 1/4 1i/2];
% approximate zero crossings:
Wo = [0.2826 0.5 0.8099 0.9073]*2*pi;
N = 4; Lth = length(Wo);      % N is the order of the allpass
W = (0:Nfft-1)'/Nfft;
Sigma = .5 + .5*cos(2*pi*W) + sin(4*pi*W);

% unweighted design
[b,a] = SingularValueAllpassSwitchCompact(Wo,N);
% weighted design
[b2,a2] = SingularValueAllpassSwitchCompact(Wo,N,sigma);

figure(1);
h = impz(b,a,1024);  h1 = circshift(h,-N+Lth/2);
h = impz(b2,a2,1024); h2 = circshift(h,-N+Lth/2);
H1 = fft(h1,Nfft); H2 = fft(h2,Nfft);
Omega = (0:Nfft-1)'/Nfft*2*pi;
plot((0:(Nfft-1))/Nfft,Sigma.*cos(unwrap(angle(H1)))); 
hold on;
plot((0:(Nfft-1))/Nfft,Sigma.*cos(unwrap(angle(H2))),'r--'); 
%plot((0:Nfft)/Nfft,(Phi+(N-Lth/2)*Omega)/pi,'r--');
grid on;
xlabel('norm. freq.');
ylabel('phase/ \pi');
legend('unweighted','weighted');
%print -depsc AllpassDesignTest5_N50.eps

disp(sprintf('integral, abs of sigma:    %f',sum(abs(Sigma))));
disp(sprintf('integral, switch:          %f',...
      sum(Sigma.*cos(unwrap(angle(H1))))));
disp(sprintf('integral, weighted switch: %f',...
      sum(Sigma.*cos(unwrap(angle(H2)))))); 
 

