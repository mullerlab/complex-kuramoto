function theta = simulate_KM( w, omega, epsilon, theta0, time, dt, method, phi )
%
% INPUT
% w - adjacency matrix (NxN)
% omega - frequencies (Nx1) (rad/s)
% epsilon - coupling strength
% theta0 - initial condition (rad)
% time - time axis (s)
% dt - timesetp (s)
% method - integration method
% phi - phase-lag
%

N = size(w,1);
theta = zeros( length(time), N ); theta(1,:) = theta0;

% integrate equations
if strcmp( method, 'euler' )
    
    for ii = 2:length(time)
        
        previous_state = theta(ii-1,:);
        
        for jj = 1:N
            dth = omega(jj) + epsilon * nansum( w(jj,:) .* sin( previous_state - previous_state(jj) - phi ), 2 );
            theta(ii,jj) = previous_state(jj) + (dth * dt);
        end
        
    end
    
elseif strcmp( method, 'ode45' )
    
    opts = odeset( 'reltol', 1e-10, 'abstol', 1e-10 );
    [ ~ , theta ] = ode45( @(t,y) KM( time, y, N, omega, epsilon, w, phi ), time, theta0, opts );
    
elseif strcmp( method, 'ode113' )
    
    opts = odeset( 'reltol', 1e-10, 'abstol', 1e-10 );    
    [ ~ , theta ] = ode113( @(t,y) KM( time, y, N, omega, epsilon, w, phi ), time, theta0, opts );
    
end

% wrap theta into [-pi, pi]
theta = angle( exp( 1i*theta ) );
