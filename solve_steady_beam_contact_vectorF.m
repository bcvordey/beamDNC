function solve_steady_beam_contact_vectorF()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) USER-PROVIDED PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M         = 500;              % Number of intervals (=> M+1 points)
    alpha     = 1;            % Beam stiffness/scale parameter
    kappa     = 2;           % Obstacle stiffness
    f_values  = [1.0,0.0,-1.0, -1.6,-2.0, -2.5, -3.0, -3.5, -4.0, -5.0]; % Vector of force values to test
    y_minus   = -0.2;            % Obstacle position
    dx        = 1.0 / M;        % Grid size
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) DISCRETIZE THE DOMAIN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xvals = linspace(0, 1, M+1)';  % Column vector of grid points
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3) LOOP OVER FORCE VALUES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_solutions = zeros(M+1, length(f_values));  % Store solutions for plotting
    
    for idx = 1:length(f_values)
        f_current = f_values(idx);
        
        % SOLVE CASE A: "No contact" assumption
        [bar_u_A, successA] = solveSystem(M, dx, alpha, kappa, f_current, y_minus, 'nocontact');
        if bar_u_A(end) > y_minus
            successA = false;
        end
        
        % SOLVE CASE B: "Contact" assumption
        [bar_u_B, successB] = solveSystem(M, dx, alpha, kappa, f_current, y_minus, 'contact');
        if bar_u_B(end) <= y_minus
            successB = false;
        end
        
        % PICK THE CONSISTENT SOLUTION
        if successA
            disp(['No-contact scenario is consistent for f = ', num2str(f_current), ...
                  ',   bar_u(1) = ', num2str(bar_u_A(end))]);
            all_solutions(:, idx) = bar_u_A;
        elseif successB
            disp(['Contact scenario is consistent for f = ', num2str(f_current), ...
                  ',   bar_u(1) = ', num2str(bar_u_B(end))]);
            all_solutions(:, idx) = bar_u_B;
        else
            warning(['No scenario is consistent for f = ', num2str(f_current), ...
                     '. Taking no-contact result as fallback.']);
            all_solutions(:, idx) = bar_u_A;  % fallback
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) PLOT THE RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure; hold on;
    cm = lines(length(f_values));  % colormap
    yline(y_minus, 'r--', 'LineWidth',1.2, 'DisplayName','Obstacle y_{-}');
    for idx = 1:length(f_values)
        plot(xvals, all_solutions(:, idx), '-', 'Color', cm(idx,:), ...
             'DisplayName',['f = ', num2str(f_values(idx))]);
    end
    
    xlabel('x');  ylabel('u(x)');
    title('Numerical solution for the Steady-State problem with f=const');
    legend('Location','southwest');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: SOLVE ONE SCENARIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bar_u, success] = solveSystem(M, dx, alpha, kappa, f_force, y_minus, mode)
    % mode: 'nocontact' or 'contact'
    % bar_u: discrete solution vector
    % success: logical indicating if assumption is consistent
    
    N   = M+1;  % number of unknowns
    A   = zeros(N, N);
    rhs = zeros(N, 1);
    
    % Indexing: bar_u(1) => \bar{u}_0, bar_u(N) => \bar{u}_M
    
    % 1) BC at x=0 => \bar{u}(0) = 0
    A(1,1) = 1;       
    rhs(1) = 0;
    
    % 2) \bar{u}_x(0) = 0 => (-3u0 + 4u1 - u2)=0 => index shift
    A(2,1) = -3;
    A(2,2) = 4;
    A(2,3) = -1;
    rhs(2) = 0;
    
    % 3) Interior points => \alpha^2 (u_{i+2} -4u_{i+1}+6u_i-4u_{i-1}+u_{i-2})/dx^4 = f_force
    for eqnIndex = 3:(N-2)
        A(eqnIndex, eqnIndex-2) =  alpha^2/(dx^4);
        A(eqnIndex, eqnIndex-1) = -4*alpha^2/(dx^4);
        A(eqnIndex, eqnIndex  ) =  6*alpha^2/(dx^4);
        A(eqnIndex, eqnIndex+1) = -4*alpha^2/(dx^4);
        A(eqnIndex, eqnIndex+2) =  alpha^2/(dx^4);
        rhs(eqnIndex) = f_force;
    end
    
    % 4) \bar{u}_{xx}(1) = 0 => u_{M-2} - 2u_{M-1} + u_M=0
    eqnIndex = N-1;
    A(eqnIndex, N-2) =  1;
    A(eqnIndex, N-1) = -2;
    A(eqnIndex, N  ) =  1;
    rhs(eqnIndex) = 0;
    
    % 5) Contact vs. no-contact BC => last row
    eqnIndex = N;
    if strcmp(mode,'nocontact')
        % no contact => \bar{u}_{xxx}(1)=0 => u_M-3u_{M-1}+3u_{M-2}-u_{M-3}=0
        A(eqnIndex, N  )   =  1/(dx^3);
        A(eqnIndex, N-1)   = -3/(dx^3);
        A(eqnIndex, N-2)   =  3/(dx^3);
        A(eqnIndex, N-3)   = -1/(dx^3);
        rhs(eqnIndex)      =  0;
    else
        % contact => - alpha^2 * bar_u_{xxx}(1)= kappa*(bar_u(M)-y_minus)
        % => alpha^2/dx^3 [u_M-3u_{M-1}+3u_{M-2}-u_{M-3}] = - kappa( u_M - y_minus )
        A(eqnIndex, N  )   =  alpha^2/(dx^3);
        A(eqnIndex, N-1)   = -3*alpha^2/(dx^3);
        A(eqnIndex, N-2)   =  3*alpha^2/(dx^3);
        A(eqnIndex, N-3)   = -1*alpha^2/(dx^3);
        % Move kappa*u_M to LHS:
        A(eqnIndex, N  )   = A(eqnIndex, N) + kappa; 
        rhs(eqnIndex)      = kappa*y_minus;
    end
    
    % Solve the system
    bar_u = A \ rhs;
    
    % We do not fully check "success" here; main script checks consistency.
    success = true;
end
