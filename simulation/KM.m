function [o] = KM( t, y, N, omega, k, w, alpha )
%
% Standard Kuramoto Model (KM)
%
% INPUT
% N         number of oscillators
% omega     vector containing instantaneous angular frequencies
% k         coupling strength
% w         N x N weight matrix (double)
%
% OUTPUT
% o         ODE output
%

o = zeros( N, 1 );
theta = y; % differential variable

for ii = 1:N
	idx = ~isnan(w(ii,:));
	o(ii) = omega(ii) + k * sum( w(ii,idx) * ( sin( theta(idx) - theta(ii) - alpha ) ) );
end
