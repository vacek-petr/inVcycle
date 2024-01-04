function approx = gs_sym(A,b,approx,numberOfIterations)
% GS_SYM Implementation of the symmetric Gauss-Seidel algorithm
%
% Inputs:
% A - matrix
% b - right-hand side vector
% approx - intial approximation of A^(-1)*b
% numberOfIterations - number of iterations to be performed
%
% Outputs:
% approx - computed approximation to A^(-1)*b

U = triu(A);
L = tril(A);
D = diag(diag(A));

for i = 1:numberOfIterations
    approx = approx + U\(D*(L\(b - A*approx)));
end
end