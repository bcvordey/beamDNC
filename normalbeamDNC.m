% 1D Beam Equation with Time Component (Finite Difference Method)
% Author: Bright Vordey
clc; clear; close all;

 % Parameters
L = 1;   % Length of the domain
T = 5;  % Total time
M = 20;   % Number of spatial points

EI = 1e-4;              % Flexural rigidity (N·m^2)
rho = 1;               % Density of the beam (kg/m)
A = 1;                 % Cross-sectional area of the beam (m^2)


% Other special parameters
beta = 0.1;
kappa = 0.1;

% Obstacle Position
obstaclePosition = -0.01;
% Velocity 
initialVelocity = 0.01;

t0 = -obstaclePosition / initialVelocity;      % Time step size

% RHS
force = -0.001;         % External forcing function f(x, t)

alpha = EI / (rho * A); % Wave propagation constant

% Spatial and time steps
dx = L / (M-1);
dt = 0.01;
N = T / dt;          % Number of time steps

% Stability condition
CFL = alpha * dt / dx^4; % Courant number (stability condition)

if CFL > 1
    error('CFL condition violated! Reduce dt or increase dx.');
end

% Tolerance
Tol = 1e-5;

% Initialize solution arrays
u = zeros(N,M);

% Initial condition
u(1,:) = 0;  
u(2,:) = 0;

% Time-stepping loop
for j = 2:N-1
    for i = 3:M-2
        u(j+1,i) = 2 * u(j,i) ...
                   - u(j-1,i) ...
                   + alpha^2 * (dt^2 / dx^4) * (u(j,i+2) - 4 * u(j,i+1) + 6 * u(j,i) - 4 * u(j,i-1) + u(j,i-2)) ...
                   + (dt^2) * force;
    end
    % Left end boundary condition
    u(j+1,1) = 0;
    u(j+1,2) = -initialVelocity * t0;


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

% Initialize grid
x = linspace(0, L, M);     % Spatial grid
t = linspace(0, T, N);     % Time grid


% Animate displacement
figure;
for n = 1:5:N
    plot(x, u(n, :), 'b', 'LineWidth', 1);
    axis([0 L -0.04 0.04]);
    xlabel('x (m)');
    ylabel('Displacement u(x,t)');
    title(['Beam Displacement at t = ' num2str((n-1)*dt) ' s']);
    grid on;
    pause(0.01);
end

% % Surface plot of displacement over time
% figure;
% surf(t, x, u', 'EdgeColor', 'none');
% xlabel('Time (s)');
% ylabel('x (m)');
% zlabel('Displacement u(x,t)');
% title('Beam Displacement Over Time');
% colormap jet;
% colorbar;
% view(45, 30);