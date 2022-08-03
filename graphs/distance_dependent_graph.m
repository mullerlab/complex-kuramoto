function A = distance_dependent_graph(N, alpha)
%
% DISTANCE-DEPENDENT GRAPH
%
% INPUT
% N - number of nodes
% alpha - power-law exponent
%
% OUTPUT
% A - adjacency matrix (weighted)
%

A = zeros( N, N );
d = nan( N );
for ii = 1:N, for jj = 1:N
        if (ii==jj), continue; end; dist = abs( ii - jj ); d(ii,jj) = min( dist, N - dist );
end; end
eta = nansum( 1 ./ d(1,:).^alpha ); A = (1.0/eta) * (1 ./ d.^alpha);
for ii = 1:N, A(ii,ii) = 0; end