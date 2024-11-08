 % Parameters
L = 1;   % Length of the domain
T = 10;  % Total time
M = 100;   % Number of spatial points
N = 10000;  % Number of time steps

% Other special parameters
beta = 0.1;
obstaclePosition = -0.5;
kappa = 0.1;

% Velocity and RHS
initialVelocity = 0.005;

force = -0.001;
t0 = -obstaclePosition / initialVelocity;
alpha = 0.0001;

% Tolerance
Tol = 1e-5;

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

    % The point before the contact point M-1;
    u(j+1,M-1)  = dt^2*(force + (2*u(j,M-1) - u(j-1,M-1))/dt^2 - (alpha^2*(24*u(j,M-1) + 20*u(j,M) - 12*u(j,M-1)))/dx^4);

    % Contact conditions at the right end
    if abs(u(j,M) - obstaclePosition) < Tol || u(j,M) < obstaclePosition

        u_2 = u(j+1,M-2);
        u_1 = u(j+1,M-1);
        u_p = u(j,M);

        contactConditionFun = @(u) alpha^2 * (((30*alpha^2*dt*u + 16*beta*dx^3*u_p - 16*beta*dx^3*u - 32*alpha^2*u_1*dt + 9*alpha^2*u_2*dt + 16*dx^3*kappa*dt*u - 16*dx^3*kappa*dt*obstaclePosition)/(7*alpha^2*dt)) ...
                                       - ((15*alpha^2*dt*u + beta*dx^3*u_p - beta*dx^3*u - 9*alpha^2*u_1*dt + alpha^2*u_2*dt + dx^3*kappa*dt*u - dx^3*kappa*dt*obstaclePosition)/(7*alpha^2*dt)) ...  
                                       + 2*u_1 - u_2) - 2*dx^3 * (kappa*u - kappa*obstaclePosition - beta*((u - u_p)/dt));

        initial_guess = u(j,M);
        u(j+1,M) = fsolve(contactConditionFun, initial_guess);  

    else
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