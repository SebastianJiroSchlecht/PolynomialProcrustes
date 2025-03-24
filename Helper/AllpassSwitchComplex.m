function [b,a] = AllpassSwitchComplex(Omegas,N);
% [b,a] = AllpassSwitchComplex(Omegas,N);
% 
%   AllpassSwitchComplex() designs a non-causal allpass filter that approximates a 
%   switching function with a gain of -1 between two specified frequencies, and a 
%   gain +1 elsewhere.
%
%   a = AllpassSwitch(Omegas,N) returns the allpass denominator coefficients in a, 
%   with Omegas = [W1 W2], with 0<= W1 < W2 <2pi defining the range over which the
%   filter gain should be negative. The parameter N is he order of the allpass. 
%   This allpass will include a delay of (N-1) samples, or needs to be advanced by 
%   (N-1) in order to have the desired phase by pi.
%
%   Input parameters:
%      Omegas [W1 W2] defines the interval where the gain should be switched
%      N      order of allpass 
%
%   Output parameter:
%      a     denominator coefficients of allpass

%  based on code by S.J. Schlecht from 14/12/2023
%  S. Weiss, UoS, 27/12/2023

%-----------------------------------------------------------------------------
% parameters   
%-----------------------------------------------------------------------------
B = Omegas(2)-Omegas(1);     % bandwidth
OmC = sum(Omegas)/2;           % centre frequency
 
%-----------------------------------------------------------------------------
% allpass design
%-----------------------------------------------------------------------------
a = AllpassSwitchReal(B,N);  % real-valued allpass centred at pi
Omega_m = OmC-pi;            % modulation
a = a.*exp(sqrt(-1)*Omega_m*(0:N)');
b = flipud(conj(a))*exp(sqrt(-1)*Omega_m);
