%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     %
% COMPLEX-VALUED APPROACH TO KURAMOTO %
% MODEL                               %
%                                     %
% BUDZINSKI ET AL. 2022               %
%                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;

% parameters
dt = 0.001; T = 4.0; t = 0:dt:T;    %time
N = 50;                             %number of oscillators
k = 25;                             %degree (2k) (ring graphs)
epsilon = 0.1;                      %coupling strength
phi = 0.0;                          %phase-lag
f_mu = 10;                          %(Hz) natural frequency
omega = ( f_mu*ones(N,1) )*2*pi;    %natural frequency

% adjacency matrix
a = ring_graph( N, k );                  %ring graph

% matrix representing the complex system
K = (epsilon) .* exp(-1i*phi) .* a;

% initial condition
theta0 = 2*pi*( rand(N,1) - 0.5 );                          %random initial conditions

% numerical simulation original KM
method = 'euler';
theta_original_km = simulate_KM( a, omega, epsilon, theta0, t, dt, method, phi );

% eigensystem
[v,d] = circulant_eigensystem( K ); %analytical eigensystem

% evaluate analytical expression
x = zeros( length(t), N ); mu = zeros( length(t), N );
win = 0.001; nwin = floor( T ./ win ); time = 0:dt:win;  %windowed approach parameters
theta_initial = theta0;                                  %initial condition for the first window
for kk = 1:nwin
    
    tmp_x = zeros( length(time), N ); tmp_mu = zeros( length(time), N );
    
    for jj = 1:length(time)
        
        %evaluation of the complex-valued model (expokit evaluation)
        tmp_x(jj,:) = exp( 1i * omega * time(jj) ) .* expv( time(jj), double(K), exp( 1i * theta_initial ) ) ;
        
        for ll = 1:N
            tmp_mu(jj,ll) = tmp_x(jj,:) * v(ll,:)'; %eigenmodes projection
        end
        
    end
    
    idx = (1:length(time)) + ((length(time)-1)*(kk-1));             %indexing time
    x(idx,:) = tmp_x;                                               %saving the solution for the complex-valued model
    theta_initial = angle( tmp_x(end,:) )';  %initial condition for the next window
    mu(idx,:) = tmp_mu;                                             %saving eigenmodes projection
    
end

theta_analytical = angle( x ); %argument of x

%fig - spatiotemporal original KM
fg1 = figure;
imagesc( 1:N, t, theta_original_km );
xlabel( 'nodes' ); ylabel( 'time (s)' ); title( 'original KM' );
set( gca, 'fontname', 'arial', 'fontsize', 18, 'linewidth', 2 );
colormap bone

%fig - spatiotemporal complex-valued model
fg2 = figure;
imagesc( 1:N, t, theta_analytical );
xlabel( 'nodes' ); ylabel( 'time (s)' ); title( 'analytical' );
set( gca, 'fontname', 'arial', 'fontsize', 18, 'linewidth', 2 );
colormap bone

% fig - eigenmodes contribution
fg3 = figure;
imagesc(1:N, t, log10(abs(mu)));
xlabel( 'modes' ); ylabel( 'time (s)' ); xlim([0.5 N+0.5]); ylim([0 T]);
set( gca, 'fontname', 'arial', 'fontsize', 18, 'linewidth', 2 );
ax = gca; ax.YDir = 'reverse';
cb = colorbar();
cb.Label.String = 'log{|\mu|}';

% fig - order parameter
fg4 = figure; hold on;
h1 = plot( t, order_parameter(theta_original_km), '-k', 'linewidth', 6 );
h2 = plot( t, order_parameter(theta_analytical), '--r', 'linewidth', 4);
le = legend( [h1 h2], {'original KM', 'analytical'} );
set( gca, 'fontname', 'arial', 'fontsize', 18, 'linewidth', 2 )
xlabel( 'time (s)' ); ylabel( '$R(t)$', 'interpreter','latex' ); ylim( [-.05 1.05] ); xlim([0 T])
