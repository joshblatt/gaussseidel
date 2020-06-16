clear; clc;

% Coeffcient matrix for Circuit 1
a1 = [0.35 -0.05 0 -0.2;
     -0.05 0.55 -0.5 0;
     0 -0.5 0.7 -0.2;
     -0.2 0 -0.2 0.44];
% Constant vector for Circuit 1
b1 = [20; 0; 0; 0];

GaussSeidelSOR(a1, b1, 0.8, 100, 0.1, 1);


% Coefficient matrix for Circuit 2
a2 = [17/105 -1/30 0 -1/10;
    -1/30 79/210 -1/5 -1/7;
    0 -1/5 23/90 -1/18;
    -1/10 -1/7 -1/18 157/315];
% Constant vector for Circuit 2
b2 = [2/7; 0; 0; 40];

GaussSeidelSOR(a2, b2, 1.5, 100, 0.1, 2);

% Coefficient matrix for Circuit 3
a3 = [11/30 -1/10;
    -1/10 19/100];
% Constant vector for Circuit 3
b3 = [16; 2.5];

GaussSeidelSOR(a3, b3, 0.8, 100, 0.1, 3);

% Prints the names of the different Voltages (ie V2, V3, etc) given the constant vector
function printVoltageNames (b)
    n = length(b);
    for i = 2: n + 1
        fprintf('      V%d  ', i)
    end
    fprintf('\n')
end

% Prints the names of the approximate errors (ie V2, V3, etc) given the constant vector
function printApproxErrorNames (b)
    n = length(b);
    for i = 2: n + 1
        fprintf('     Ea%d  ', i)
    end
    fprintf('\n')
end

% Calculates the currents given the voltages and circuit number and prints out the result
function calculateCurrents(x, circuitNum)
    % Create variables for indicies of each voltage (to make calculations clearer)
    V2 = 1;
    V3 = 2;
    V4 = 3;
    V5 = 4;
    
    if circuitNum == 1
        % currents = [i12 i52 i32 i65 i54 i43]
        currents = [0 0 0 0 0 0];
        % Use Ohm's law to calculate currents based on given circuit
        currents(1) = (200 - x(V2)) / 10;
        currents(2) = (x(V5) - x(V2)) / 5;
        currents(3) = (x(V3) - x(V2)) / 20;
        currents(4) = (-x(V5)) / 25;
        currents(5) = (x(V5) - x(V4)) / 5;
        currents(6) = (x(V4) - x(V3)) / 2;
        fprintf('     i12       i52       i32       i65       i54       i43\n');
        disp(currents);
    elseif circuitNum == 2
        % currents = [i12 i52 i32 i65 i54 i43 i53]
        currents = [0 0 0 0 0 0 0];
        % Use Ohm's law to calculate currents based on given circuit
        currents(1) = (10 - x(V2)) / 35;
        currents(2) = (x(V5) - x(V2)) / 10;
        currents(3) = (x(V3) - x(V2)) / 30;
        currents(4) = (200 - x(V5)) / 5;
        currents(5) = (x(V5) - x(V4)) / 18;
        currents(6) = (x(V4) - x(V3)) / 5;
        currents(7) = (x(V5) - x(V3)) / 7;
        fprintf('     i12       i52       i32       i65       i54       i43       i53\n');
        disp(currents);
    elseif circuitNum == 3
        % currents = [i12 i25 i23 i35 i43]
        currents = [0 0 0 0 0];
        % Use Ohm's law to calculate currents based on given circuit
        currents(1) = (80 - x(V2)) / 5;
        currents(2) = x(V2) / 15;
        currents(3) = (x(V2) - x(V3)) / 10;
        currents(4) = x(V3) / 25;
        currents(5) = (50 - x(V3)) / 20;
        fprintf('     i12       i25       i23       i35       i43\n');
        disp(currents);
    end
end

% Calculates voltage vector using linsolve, calculates currents, and prints result
function linsolveEquivalent(A, b, circuitNum)
    fprintf('MATLAB Equation Solver Equivalent\n');
    x = linsolve(A, b);
    printVoltageNames(x);
    disp(x');
    calculateCurrents(x, circuitNum);
end

% Function that does Gauss Seidel SOR given A (coefficient matrix)
% B (constant vector), lambda (where 0 < lambda < 2), maxIterations,
% maxApproximateError (%), and the circuitNum
function GaussSeidelSOR(A, b, lambda, maxIterations, maxApproxError, circuitNum)
    n = length(b); % Get length of coefficient matrix, save in n
    % x is the unknown solution vector (represented as an nxn matrix)
    x = zeros(n);
    % x_0 is the previous iteration of x (used for approx error calculations)
    x_0 = zeros(n);
    
    % D is a matrix with the diagonal elements from A, with all other elements being 0
    D = diag(diag(A));
    % Strictly upper triangular matrix of A
    U = triu(A,1);
    % Strictly lower triangular matrix of A
    L = tril(A, -1);
    % d is the inverse of the diagonal matrix + lambda * lower triangular
    % (used in SOR)
    d = inv(D + lambda*L);
    
    % Iterate from 1 - maxIterations (at most)
    for i = 1: maxIterations
        % Find the newest iteration of x based on L9.2 slide 6
        x = d * (lambda * b - (lambda * U + D * (lambda - 1)) * x_0);
        % Calculate approximate error
        ea = abs((diag(x) - diag(x_0)) * 100./ diag(x));
        % If the error is less than the max approximate error
        if max(ea) < maxApproxError
            
            % Exit the loop (answer has been found)
            break;
        end
        % Set previous value of x to x for next iteration
        x_0 = x;
    end
    
    fprintf('Circuit: %d\n', circuitNum);
    fprintf('Iteration Number %d \n', i);
    fprintf('Lambda: %.1f\n', lambda);
    
    % Prints voltages + labels
    printVoltageNames(x);
    disp(diag(x)');
    
    % Prints approximate errors + labels
    printApproxErrorNames(x);
    disp(ea');
    
    % Calculates and prints currents
    calculateCurrents(diag(x)', circuitNum);
    
    % Calculates voltages using linsolve and prints result
    linsolveEquivalent(A, b, circuitNum);
end