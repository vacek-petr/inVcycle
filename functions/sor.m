function approx = sor(A,b,approx,numberOfIterations,dampingParameter) 
% SOR Implementation of the SOR algorithm
%
% Inputs:
% A - matrix
% b - right-hand side vector
% approx - intial approximation of A^(-1)*b
% numberOfIterations - number of iterations to be performed
% dampingParameter
%
% Outputs:
% approx - computed approximation to A^(-1)*b


D = diag(diag(A));
U = triu(A, 1);
L = tril(A, -1);

for i = 1:numberOfIterations
     approx = (D + dampingParameter*L) \ (dampingParameter*b - (dampingParameter*U + (dampingParameter-1)*D)*approx);
end
