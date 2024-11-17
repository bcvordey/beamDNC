% 1D Beam Equation with Time Component (Finite Difference Method)
% Author: Bright Vordey
clc; clear; close all;

 % Parameters
L = 1;   % Length of the domain
T = 5;  % Total time
M = 50;   % Number of spatial points


EI = 1e-4;              % Flexural rigidity (N·m^2)
rho = 1;               % Density of the beam (kg/m)
A = 1;                 % Cross-sectional area of the beam (m^2)


% Other special parameters
beta = 0.1;
obstaclePosition = 0;
kappa = 0.1;

% Velocity and RHS
initialVelocity = 1;
force = -0.001;         % External forcing function f(x, t)


t0 = -obstaclePosition / initialVelocity;
alpha = EI / (rho * A); % Wave propagation constant

% Spatial and time steps
dx = L / (M-1);
dt = 0.001;            % Time step size
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
    u(j+1,1) = u(j+1,2);
    u(j+1,M-1) = 0;            % Enforce clamped boundary
    u(j+1,M) = 0;              % u(L, t) = 0
end

% Initialize grid
x = linspace(0, L, M);     % Spatial grid
t = linspace(0, T, N);     % Time grid


% Animate displacement
figure;
for n = 1:5:N
    plot(x, u(n, :), 'b', 'LineWidth', 1);
    axis([0 L -0.1 0.1]);
    xlabel('x (m)');
    ylabel('Displacement u(x,t)');
    title(['Beam Displacement at t = ' num2str((n-1)*dt) ' s']);
    grid on;
    pause(0.01);
end

% Surface plot of displacement over time
figure;
surf(t, x, u', 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('x (m)');
zlabel('Displacement u(x,t)');
title('Beam Displacement Over Time');
colormap jet;
colorbar;
view(45, 30);