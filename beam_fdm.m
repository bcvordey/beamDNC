% Parameters
L = 1;   % Length of the domain
T = 10;  % Total time
M = 10;   % Number of spatial points
N = 10000;  % Number of time steps

% Other special parameters
beta = 0.1;
obstaclePosition = -2;
kappa = 0.01;

% Velocity and RHS
initialVelocity = 1;
force = -10;
t0 = - obstaclePosition / initialVelocity;

alpha = 0.0001;

% Tolerance
Tol = 1e-3;

% Spatial and time steps
dx = L / M;
dt = (T - t0) / N;

% Initialize solution arrays
u = zeros(N,M);

% Initial condition
u(1,:) = 0;  
u(2,:) = initialVelocity * dt;

% Time-stepping loop
for j = 2:N-1
    for i = 3:M-2
        u(j+1,i) = 2 * u(j,i) ...
                   - u(j-1,i) ...
                   + alpha^2 * (dt^2 / dx^4) * (u(j,i+2) - 4 * u(j,i+1) + 6 * u(j,i) - 4 * u(j,i-1) + u(j,i-2)) ...
                   + (dt^2) * force;
    end
    
    % Left end boundary condition
    u(j+1,1) = u(j+1,2);

    % Contact conditions at the right end
    if abs(u(j,M) - obstaclePosition) < Tol || u(j,M) < obstaclePosition

        u(j+1,M-1) = - ((dt * alpha^2) / (2 * dx^2 * beta)) * (2 * u(j,M-2) - u(j,M-3) - u(j,M)) ...
                     - ((dt * kappa) / beta) * (u(j,M-1) - obstaclePosition) ...
                     + u(j-1,M-1);

        u(j+1,M) = u(j+1,M-1);
    else
        u(j+1,M-1) = u(j+1, M-2);
        u(j+1,M) = u(j+1, M-1);
    end
end

disp(u);

% Plot results
x_values = linspace(0, L, M);  
figure;
plot(x_values, u(N,:), 'r', 'DisplayName', 'Left end');        
yline(obstaclePosition, '--g', 'DisplayName', 'Obstacle');  % Obstacle
xlabel('Time (t)');
ylabel('u');
title('FDM Solution');
legend show;
hold off;
