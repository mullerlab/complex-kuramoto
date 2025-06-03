%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
% OPERATOR EXPRESSION FOR THE KURAMOTO MODEL    %
% BUDZINSKI ET AL. (CHAOS 2022)                 %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup
clearvars; clc;

% parameters
dt = 0.001; T = 1.0; t = 0:dt:T;         %time
N = 50;                                  %number of oscillators
epsilon = 0.75;                          %coupling strength
phi = 0.0;                               %phase-lag
f_mu = 5;                                %(Hz) natural frequency
omega = ( f_mu*ones(N,1) )*2*pi;         %natural frequency

% adjacency matrix
k = 10;                                  %degree (2k)
a = ring_graph( N, k );                  %ring graph

% aggregate matrix and analytical eigenspectrum
K = epsilon .* exp(-1i*phi) .* a;
[v,d] = circulant_eigensystem( K );

% initial condition
rng(1); theta0 = 2*pi*( rand(N,1) - 0.5 ); %random initial conditions (w/RNG seed)

% numerical simulation - Kuramoto equations
theta_km = simulate_KM( a, omega, epsilon, theta0, t, dt, phi );

% evaluate operator expression
x = nan( length(t), N ); x(1,:) = exp( 1i * theta0 ); tau = 0.001;
propagate = exp( 1i * omega * tau ) .* expm( tau * double(K) ); %propagator

for jj = 2:length(t)

    x(jj,:) = propagate * reshape( x(jj-1,:), [], 1 );
    x(jj,:) = x(jj,:) ./ abs( x(jj,:) ); % unit modulus

end
theta_cv = angle( x ); % argument of x

% evaluate eigenmodes
mu = nan( length(t), N );

for jj = 1:length(t)
    
    for ii = 1:N, mu(jj,ii) = v(:,ii)' * x(jj,:).'; end

end

% fig - spatiotemporal dynamics
fg1 = figure; set( fg1, 'position', [1   688   723   310] )
subplot(121); h1 = imagesc( 1:N, t, theta_km ); colormap hsv; title( 'kuramoto' )
set( gca, 'linewidth', 2, 'fontsize', 18 ); xlabel( 'nodes' ); ylabel( 'time (s)' )
subplot(122); h2 = imagesc( 1:N, t, theta_cv ); colormap hsv; title( 'cv approach' )
set( gca, 'linewidth', 2, 'fontsize', 18 ); xlabel( 'nodes' ); ylabel( 'time (s)' )

% fig - order parameter
fg2 = figure; hold on;
set( fg2, 'position', [92   410   560   198] )
plot( t, order_parameter(theta_km), '-k', 'linewidth', 6 );
plot( t, order_parameter(theta_cv), '--r', 'linewidth', 4);
le = legend({'kuramoto', 'cv approach'} );
set( gca, 'fontname', 'arial', 'fontsize', 18, 'linewidth', 2 )
xlabel( 'time (s)' ); ylabel( '$R(t)$', 'interpreter','latex' ); ylim( [-.05 1.05] ); xlim([0 T])

% fig - eigenmodes
fg3 = figure; 
set( fg3, 'position', [114    17   495   313] )
imagesc( 1:N, t, abs(mu));
cb = colorbar();
cb.Label.String = '|\mu|';
set( gca, 'linewidth', 2, 'fontsize', 18 ); xlabel( 'nodes' ); ylabel( 'time (s)' )


function A = ring_graph( N, k )
%
% RING GRAPH
%
% INPUT
% N - number of nodes
% k - node degree
%
% OUTPUT
% A - adjacency matrix
%

A = false( N, N );
m = -k:k;m( m == 0 ) = [];

% generate ring graph
for ii = 0:(N-1)
    for jj = m

        A( ii+1, mod(ii+jj,N)+1 ) = 1;

    end
end

end

function [v,d] = circulant_eigensystem( a )
%
% CIRCULANT EIGENSYSTEM
%
% INPUT
% a - adjacency matrix (has to be circulant)
%
% OUTPUT
% v - eigenvectors (NxN matrix)
% d - eigenvalues (NxN diagonal matrix)
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

end

function theta = simulate_KM( w, omega, epsilon, theta0, time, dt, phi )
%
% NUMERICAL INTEGRATION OF THE KURAMOTO MODEL
%
% INPUT
% w - adjacency matrix (NxN)
% omega - frequencies (Nx1) (rad/s)
% epsilon - coupling strength
% theta0 - initial condition (rad)
% time - time axis (s)
% dt - timesetp (s)
% phi - phase-lag
%

N = size(w,1);
theta = zeros( length(time), N ); theta(1,:) = theta0;
    
for ii = 2:length(time)

    previous_state = theta(ii-1,:);

    for jj = 1:N
        dth = omega(jj) + epsilon * nansum( w(jj,:) .* sin( previous_state - previous_state(jj) - phi ), 2 );
        theta(ii,jj) = previous_state(jj) + (dth * dt);
    end

end

% wrap theta into [-pi, pi]
theta = angle( exp( 1i*theta ) );

end

function r = order_parameter( ang )
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

end
