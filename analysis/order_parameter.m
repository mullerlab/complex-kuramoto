function r = order_parameter( ang, NN )
%
% ORDER PARAMETER
%
% INPUT
% ang - TxN matrix of angular values
%
% OUTPUT
% r - order parameter \in [0,1]
%

NN = size(ang,2);
r = (1/NN) * sum( exp(1i*ang), 2 );
r = abs(r);
