%% 
clc;close all;clear;

% This code solves the free vibration of a 1D beam using Finite Difference 

%% User Data
L=1; % total length of the beam (m)
dx=0.025; %step size for the length (m)
x=0:dx:L; % locations of nodes
M=length(x); % number of nodes


animationSpeed = 1e4;
timeIntervalAdjust = 10;



T=2; % final computation time
dt=4e-6; % time increment
t=0:dt:T; % time points
N=length(t); % number of time steps used
Tol = 1e-10;

% % Special Parameters
beta = 0.01;
kappa = 1;
obstacle = -0.02;
force = -1.65;

%Other parameters
%beam with rectangular cross section
A=5*2; % m^2
rho=2700; % kg/m
E=85e9; % Pa . %%%% 100 for beta 0.01 and 65 for beta 0.1
I=0.2*(0.05^3)/12; % m^4


alpha=((rho*A)/(E*I))^0.5;

disp('lambda must be less than or equal 1 for stability')
lambda = dt/((dx^2)*(alpha));

tau1 = alpha^2/dx^3;
tau2 = beta/dt;

phi1 = 2*tau1 + kappa - tau2;
phi2 = -tau1 - 2*kappa + 2*tau2;



%% Boundary Conditions and Initial Condition

u=zeros(M,N); % intialization of displacements for all nodes at all time points

% Boundary conditions
u(1,:)=0;  
u(2,:)=u(1,:); 


% initial conditions
P=3*E*I*0.1/(L^3); % force for a 10-cm upward displacement at the end
u(:,1)=P*(x.^2).*(3*L-x)/(6*E*I); % initial displacements at all nodes, q(x)
u(:,2)=u(:,1);

ymax=max(abs(u(:,1)));

%% Solution

for j=3:N   % loop over time points
    for i=3:M-2  % loop over nodes
        u(i,j)= -(lambda^2)* ...
                 (u(i-2,j-1)-4*u(i-1,j-1)-4*u(i+1,j-1)+u(i+2,j-1)+6*u(i,j-1)) ...
                 -u(i,j-2)+2*u(i,j-1) + (dt^2) * force;        
    end

    % Boundary conditions
   
    if abs(u(M,j-1) - obstacle) < Tol || u(M,j-1) < obstacle

        u(M-1,j)= (tau1 * u(M-3,j) - phi1 * u(M-2,j) + tau2 * u(M,j-1) - kappa * obstacle)/(phi2);

        u(M,j)= 2 * u(M-1,j) - u(M-2,j);

    else
        u(M-1,j)=2*u(M-2,j)-u(M-3,j);
        u(M,j)=3*u(M-2,j)-2*u(M-3,j);
    end
       
end    


%% visualization
h=plot(x,u(:,1),'-','LineWidth',5);
xlabel('x');
ylabel('displacements u(x,t)');
title('vibrations of a beam with Obstacle')
axis([0 1.1*L -1.1*ymax 1.1*ymax]);

% Add obstacle at x = 1 with height from -0.1 to -0.05
hold on; % Allow additional graphics on the same axes
line([1, 1], [-1.1*ymax, obstacle], 'Color', 'k', 'LineWidth', 5); % Vertical bar
grid on;
hold off;

for kk=2:animationSpeed:N
    set(h,'YData',u(:,kk));
    drawnow; 
end

% --- Sub-Sample Time and Data ---
idx        = 1:timeIntervalAdjust:length(t);
t_sub      = t(idx);
dist0M     = u(M,:);
dist0M_sub = dist0M(idx);

% --- Create Figure & Plot ---
figure;
plot(t_sub, dist0M_sub, ...
     '-', ...                    % a line with markers at each sub-sample
     'LineWidth', 0.5, ...        % thicker line for clarity
     'MarkerSize', 1, ...         % small markers
     'DisplayName', 'Beam Position');
hold on;                          % so we can add more plot elements easily

% Add a horizontal line for the obstacle
yline(obstacle, 'r--', ...
      'DisplayName', 'Obstacle Position');

% --- Axis Labels & Title ---
xlabel('Time (s)');
ylabel('u(M,t)');

% You can split the title into multiple lines for readability
titleStr = sprintf([ ...
    'Position of the Contacting End of the Beam\n', ...
    '\\alpha^2=%.2f, \\beta=%.2f, \\kappa=%.3f, obstacle=%.2f, force=%.2f' ...
    ], alpha^2, beta, kappa, obstacle, force);
title(titleStr);

% --- Cosmetics & Legend ---
grid on;
legend('Location','northeast');
hold off;