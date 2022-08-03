function [v,d] = circulant_eigensystem( a )
%
% CIRCULANT EIGENSYSTEM
%
% INPUT
% a - adjacency matrix (YOU must verify circulant)
%
% OUTPUT
% v - eigenvectors (NxN matrix)
% d - eigenvalues (NxX diagonal matrix)
%

N = size( a, 1 );

v = zeros( N ); d = zeros( N );
for ii = 1:N
    for jj = 1:N
        v(ii,jj) = (1/sqrt(N)) * exp( -2*pi*1i/N * (ii-1) * (jj-1) );
    end
end

for ii = 1:N
    d(ii,ii) = sum( a(1,:) .* exp( -2*pi*1i/N * ( ii - 1 ) * ( (1:N) - 1 ) ) );
end
