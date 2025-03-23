function singlesolve_steady_beam_contact()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) USER PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M       = 200;       % Number of intervals (so M+1 points)
    alpha   = 1;      % Positive constant (scaled)
    kappa   = 1.0;     % Contact stiffness
    f       = -1.6;      % Constant force
    y_minus = -0.2;      % Obstacle position
    dx      = 1.0 / M;  % Grid size

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) DISCRETIZATION SETUP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We'll create an initial guess or set of unknowns:
    %   bar_u(i) ~ bar(u) at x_i, i=0,...,M
    % The domain is x in [0,1].
    xvals = linspace(0, 1, M+1)';  % Column vector

    % We will handle boundary conditions:
    %   bar_u(0)=0
    %   bar_u_x(0)=0
    %   bar_u_xx(1)=0
    %   -alpha^2 bar_u_xxx(1)= kappa*(bar_u(1)-y_minus)_-

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3) SOLUTION STRATEGY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We'll do a 2-case approach:
    % CASE A: No contact => (bar_u(M) >= y_minus)
    %   Then (bar_u(M)-y_minus)_- = 0 => bar_u(M) >= y_minus
    %   => bar_u_xxx(1)=0
    %
    % CASE B: Contact => (bar_u(M) < y_minus)
    %   => -alpha^2 bar_u_xxx(1) = kappa*(bar_u(M)-y_minus)
    %
    % We'll guess a case, solve the linear system, then check consistency.
    % If inconsistent, we switch the case and re-solve.

    % Start with an initial guess (all zeros or linear).

    % We solve iteratively: try case A, check bar_u(M), if conflict => case B.
    % Then solve case B, check bar_u(M). Possibly no iteration needed if one pass is consistent.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) SOLVE CASE A: NO CONTACT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [bar_u_A, successA] = solveSystem(M, dx, alpha, kappa, f, y_minus, 'nocontact');

    % Check consistency => if bar_u_A(M+1) < y_minus => case A fails
    if bar_u_A(M+1) < y_minus
        successA = false;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5) SOLVE CASE B: CONTACT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [bar_u_B, successB] = solveSystem(M, dx, alpha, kappa, f, y_minus, 'contact');

    % Check consistency => if bar_u_B(M+1) >= y_minus => case B fails
    if bar_u_B(M+1) >= y_minus
        successB = false;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6) PICK SOLUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if successA
        disp('No-contact case is consistent with bar_u(1) >= y_minus.');
        bar_u_sol = bar_u_A;
    elseif successB
        disp('Contact case is consistent with bar_u(1) < y_minus.');
        bar_u_sol = bar_u_B;
    else
        % In more complicated scenarios, you might refine or iterate further,
        % but typically one of these cases should hold.
        warning('Neither contact nor no-contact is consistent. Check parameters or refine approach.');
        % For demonstration, pick whichever is physically more plausible.
        bar_u_sol = bar_u_A;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 7) PLOT THE SOLUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;  hold on;
    plot(xvals, bar_u_sol, '-','LineWidth',1);
    yline(y_minus,'r--','LineWidth',1.2);
    xlabel('x');
    ylabel('bar\_u(x)');
    legend('Steady-State Displacement','Obstacle Position','Location','Best');
    title('Steady-State Beam Displacement with Contact Condition');
    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Function: solveSystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bar_u, success] = solveSystem(M, dx, alpha, kappa, f, y_minus, mode)
    % mode is 'nocontact' or 'contact'
    % Returns bar_u if the system solves, success= true/false if solution is feasible

    % Number of unknowns = M+1
    N = M+1;

    % We'll build a matrix A and right-hand side rhs for A * bar_u = rhs
    A   = zeros(N, N);
    rhs = zeros(N, 1);

    % Indices: i=1..N => corresponds to x_{i-1} in [0,M]
    %   bar_u(1) ~ bar_u_0, bar_u(N) ~ bar_u_M

    % 1) BC: bar_u(0)=0 => bar_u(1) = 0
    A(1,1) = 1;   % bar_u(0)=0
    rhs(1)  = 0;

    % 2) BC: bar_u_x(0)=0 => derivative=0 => central difference approx
    %   -3u0 + 4u1 - u2 = 0
    % index shift: i=2 => corresponds to eqn for slope at x=0
    A(2,1) = -3;
    A(2,2) = 4;
    A(2,3) = -1;
    rhs(2) = 0;

    % 3) Interior eqns: alpha^2 (u_{i+2}-4u_{i+1}+6u_i -4u_{i-1}+u_{i-2})/dx^4 = f
    % i ranges from 3..(N-2)
    % note we have to carefully shift by 2 indices
    for eqnIndex = 3 : (N-2)
        i = eqnIndex;  % i in matrix sense
        % alpha^2/dx^4 * [u_{i+2}-4u_{i+1}+6u_i -4u_{i-1}+u_{i-2}]
        A(eqnIndex, i-2) = alpha^2/(dx^4);
        A(eqnIndex, i-1) = -4*alpha^2/(dx^4);
        A(eqnIndex, i   ) = 6*alpha^2/(dx^4);
        A(eqnIndex, i+1 ) = -4*alpha^2/(dx^4);
        A(eqnIndex, i+2 ) = alpha^2/(dx^4);
        rhs(eqnIndex) = f;  % constant
    end

    % 4) BC: bar_u_xx(1)=0 => [u_{M-2} - 2u_{M-1} + u_M] = 0
    % This is eqn for eqnIndex=(N-1)
    eqnIndex = N-1;
    A(eqnIndex, N-2) = 1;   % bar_u_{M-2}
    A(eqnIndex, N-1) = -2;  % bar_u_{M-1}
    A(eqnIndex, N  ) = 1;   % bar_u_{M}
    rhs(eqnIndex) = 0;

    % 5) BC: contact or no-contact at x=1
    eqnIndex = N;
    if strcmp(mode, 'nocontact')
        % CASE A: no contact => bar_u(1)>=y_minus => (bar_u(1)-y_minus)_- = 0 => bar_u_xxx(1)=0
        % bar_u_xxx(1) ~ [u_M - 3u_{M-1} + 3u_{M-2} - u_{M-3}]/dx^3 = 0
        A(eqnIndex, N  )   =  1/(dx^3);   % bar_u_M
        A(eqnIndex, N-1)   = -3/(dx^3);   % bar_u_{M-1}
        A(eqnIndex, N-2)   =  3/(dx^3);   % bar_u_{M-2}
        A(eqnIndex, N-3)   = -1/(dx^3);   % bar_u_{M-3}
        rhs(eqnIndex)      = 0;
    else
        % CASE B: contact => bar_u(1)<y_minus => -alpha^2 bar_u_xxx(1)= kappa*(bar_u(1)-y_minus)
        % bar_u_xxx(1) ~ [u_M - 3u_{M-1} + 3u_{M-2} - u_{M-3}]/dx^3
        % => alpha^2/dx^3*[u_M-3u_{M-1}+3u_{M-2}-u_{M-3}] = -kappa*(u_M-y_minus)
        %
        % => alpha^2/dx^3*u_M -3*... = -kappa u_M + kappa*y_minus
        % We'll gather terms:

        % left side
        A(eqnIndex, N  )   =  alpha^2/(dx^3);
        A(eqnIndex, N-1)   = -3*alpha^2/(dx^3);
        A(eqnIndex, N-2)   =  3*alpha^2/(dx^3);
        A(eqnIndex, N-3)   = -1*alpha^2/(dx^3);

        % Move -kappa*u_M to left side => +kappa on diagonal of bar_u_M
        A(eqnIndex, N  )   = A(eqnIndex, N) + kappa;  % combine terms
        rhs(eqnIndex)      =  kappa*y_minus;          % constant from +kappa*y_minus
    end

    % Now we have a (N x N) system: A * bar_u = rhs
    % Solve
    bar_u = A \ rhs;

    % For success check or consistency, we do minimal checks here
    success = true;  % We'll check after return if it matches the sign condition
end
