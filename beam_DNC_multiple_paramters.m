clc; close all; clear;

% This code solves the free vibration of a 1D beam with DNC condition using a Finite Difference Method 
% and loops over multiple values of kappa, beta, and force (f).

%% User-defined vectors for parametric study
kappaVec = [0.1, 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6];
betaVec  = 0.01;
forceVec = [-5, -10];

%% Fixed user data (remains the same for each combination)
L  = 1;      % total length of the beam (m)
dx = 0.025;  % step size for the length (m)
x  = 0:dx:L; % locations of nodes
M  = length(x); % number of nodes

animationSpeed     = 1e5;  % controls how quickly the shape is updated
timeIntervalAdjust = 10;    % for subsampling time

T   = 2;       % final computation time
dt  = 4e-6;    % time increment
t   = 0:dt:T;  
N   = length(t);

Tol       = 1e-10;
obstacle  = -0.02;    % obstacle location (in displacement)

% Beam physical parameters
A   = 5*2;       % cross-sectional area (m^2)
rho = 2700;      % density (kg/m)
E   = 85e9;      % Young's modulus (Pa)
I   = 0.2*(0.05^3)/12; % second moment of area (m^4)

% Precompute alpha (depends on beam geometry and materials, not on kappa/beta/f)
alpha = sqrt((rho*A)/(E*I));

% Check stability parameter
lambda = dt / ((dx^2)*(alpha));
disp(['lambda must be <= 1 for stability. Current lambda = ', num2str(lambda)]);

%% Nested loops over kappa, beta, force
for kk = 1:length(kappaVec)
    for bb = 1:length(betaVec)
        for ff = 1:length(forceVec)

            kappa = kappaVec(kk);
            beta  = betaVec(bb);
            f     = forceVec(ff);

            %--- Additional parameters for the contact boundary condition
            tau1 = alpha^2 / dx^3;
            tau2 = beta / dt;
            phi1 = 2*tau1 + kappa - tau2;
            phi2 = -tau1 - 2*kappa + 2*tau2;

            %% Initialize displacement array for all nodes & time steps
            u = zeros(M, N);

            % Boundary conditions
            % Node 1 and 2 are always zero here
            u(1,:) = 0;
            u(2,:) = 0;

            % Initial condition: (10 cm = 0.1 m) upward displacement shape
            % "P" from beam bending formula for an end deflection of 0.1 m
            P       = 3*E*I*0.1/(L^3);     
            u(:,1)  = P*(x.^2).*(3*L - x)/(6*E*I);
            u(:,2)  = u(:,1);

            % Helpful for axis scaling later
            ymax = max(abs(u(:,1)));

            %% Main FDM Loop
            for j = 3:N
                for i = 3:M-2 
                    u(i,j) = ...
                         -(lambda^2)*( ...
                            u(i-2,j-1) - 4*u(i-1,j-1) - 4*u(i+1,j-1) + u(i+2,j-1) + 6*u(i,j-1) ...
                          ) ...
                         - u(i,j-2) ...
                         + 2*u(i,j-1) ...
                         + (dt^2)*f;
                end

                % -- Contact boundary condition at the right end (i = M)
                if abs(u(M,j-1) - obstacle) < Tol || (u(M,j-1) < obstacle)
                    % The beam tip has reached or crossed the obstacle
                    u(M-1,j) = (tau1 * u(M-3,j)  - phi1 * u(M-2,j)  + tau2 * u(M,j-1) ...
                               - kappa * obstacle ) / phi2;
                    u(M,  j) = 2*u(M-1,j) - u(M-2,j);
                else
                    % No contact
                    u(M-1,j) = 2*u(M-2,j) - u(M-3,j);
                    u(M,  j) = 3*u(M-2,j) - 2*u(M-3,j);
                end
            end

            % %% Visualization #1: Animate shape of the beam
            % figure('Name', ...
            %     sprintf('BeamShape_beta%.2g_kappa%.2g_force%.2g', beta, kappa, f), ...
            %     'Color','w');
            % 
            % h = plot(x, u(:,1), '-', 'LineWidth', 5);
            % hold on;
            % 
            % % Add the obstacle line at x = L
            % line([L, L], [-1.1*ymax, obstacle], 'Color', 'k', 'LineWidth', 5);
            % 
            % xlabel('x');
            % ylabel('Displacement u(x,t)');
            % title({
            %     'Vibration of a Beam with Obstacle', ...
            %     sprintf('\\beta=%.2g, \\kappa=%.2g, f=%.2g', beta, kappa, f)
            %     });
            % axis([0, 1.1*L, -1.1*ymax, 1.1*ymax]);
            % grid on;
            % hold off;
            % 
            % % Optional animation over time
            % for jj = 2:animationSpeed:N
            %     set(h, 'YData', u(:,jj));
            %     drawnow; 
            % end

            %% Visualization #2: Contact end displacement over time
            % Subsample time and displacement for a simpler plot
            idx        = 1:timeIntervalAdjust:length(t);
            t_sub      = t(idx);
            dist0M     = u(M,:);       % displacement at the right end
            dist0M_sub = dist0M(idx);

            fig2Name = sprintf('ext_beta%.2g_kappa%.2g_force%.2g', beta, kappa, f);
            fig2 = figure('Name', fig2Name, 'Color','w');

            plot(t_sub, dist0M_sub, '-', ...
                 'LineWidth', 0.5, 'MarkerSize', 1, ...
                 'DisplayName', 'Beam Tip Position');
            hold on;

            % Obstacle line
            yline(obstacle, 'r--', 'DisplayName', 'Obstacle Position');

            xlabel('Time (s)');
            ylabel('u(M,t)');
            title({
                'Contact End of the Beam vs. Time', ...
                sprintf('\\alpha^2=%.4f, \\beta=%.2g, \\kappa=%.2g, obstacle=%.2f, force=%.2f', ...
                    alpha^2, beta, kappa, obstacle, f)
                });
            grid on;
            legend('Location','northeast');
            hold off;

            % Save the figure
            saveas(fig2, [fig2Name, '.jpg']);

        end
    end
end
