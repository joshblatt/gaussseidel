clear; clc;

% Coeffcient matrix for Circuit 1
a1 = [0.35 -0.05 0 -0.2;
     -0.05 0.55 -0.5 0;
     0 -0.5 0.7 -0.2;
     -0.2 0 -0.2 0.44];
% Constant vector for Circuit 1
b1 = [20; 0; 0; 0];

V1 = GaussSeidelSOR(a1, b1, 0.8, 100, 0.1);
printFinalResult(V1);
disp((linsolve(a1, b1))');

% Coefficient matrix for Circuit 2
a2 = [17/105 -1/30 0 -1/10;
    -1/30 79/210 -1/5 -1/7;
    0 -1/5 23/90 -1/18;
    -1/10 -1/7 -1/18 157/315];
% Constant vector for Circuit 2
b2 = [2/7; 0; 0; 40];

V2 = GaussSeidelSOR(a2, b2, 1.5, 100, 0.1);
printFinalResult(V2);
disp((linsolve(a2, b2))');

% Coefficient matrix for Circuit 3
a3 = [11/30 -1/10;
    -1/10 19/100];
% Constant vector for Circuit 3
b3 = [16; 2.5];

V3 = GaussSeidelSOR(a3, b3, 0.8, 100, 0.1);
printFinalResult(V3);

% Prints the names of the different Voltages (ie V2, V3, etc)
% given the constant vector
function printVoltageNames (b)
    n = length(b);
    for i = 2: n + 1
        fprintf('      V%d  ', i)
    end
    fprintf('\n')
end

% Prints the names of the approximate errors (ie V2, V3, etc)
% given the constant vector
function printApproxErrorNames (b)
    n = length(b);
    for i = 2: n + 1
        fprintf('     Ea%d  ', i)
    end
    fprintf('\n')
end

% Prints the voltages and currents
function printFinalResult(Voltages)
    fprintf('Final Result\n');
    printVoltageNames(Voltages);
    disp(diag(Voltages)');
    
%     if circuitNum == 1
%         currents = 1:6;
%     end
end

% Function that does Gauss Seidel SOR given A (coefficient matrix)
% B (constant vector), lambda, maxIterations, and the maxApproximateError
function [x] = GaussSeidelSOR(A, b, lambda, maxIterations, maxApproxError)
    n = length(b); % Get length of coefficient matrix, save in n
    % x is the unknown solution vector (represented as an nxn matrix)
    x = zeros(n);
    % x_0 is the previous iteration of x (used for approx error
    % calculations)
    x_0 = zeros(n);
    
    % D is a matrix with the diagonal elements from A, with all other
    % elements being 0
    D = diag(diag(A));
    % Strictly upper triangular matrix of A
    U = triu(A,1);
    % Strictly lower triangular matrix of A
    L = tril(A, -1);
    % d is the inverse of the diagonal matrix + lambda * lower triangular
    % (used in SOR)
    d = inv(D + lambda*L);
    
    % Iterate from 1 - maxIterations
    for i = 1: maxIterations
        % Find the newest iteration of x based on L9.2 slide 6
        x = d * (lambda * b - (lambda * U + D * (lambda - 1)) * x_0);
        %fprintf('Iteration Number %d \n', i);
        %printVoltageNames(b);
        % Display x as a horizontal vector
        %disp(diag(x)');
        % Calculate approximate error
        ea = abs((diag(x) - diag(x_0)) * 100./ diag(x));
        %printApproxErrorNames(b);
        % Display approximate errors as horizontal vector
        %disp(ea')
        % If the error is less than the max approximate error
        if max(ea) < maxApproxError
            % Exit the loop (answer has been found)
            break;
        end
        % Set previous value of x to x for next iteration
        x_0 = x;
    end
end