function estimate = computeErrAnormMLEstimate(A,P,RHS,J,approx,options)
% COMPUTEERRANORMMLESTIMATE implementation of the multilevel error
% estimator on the Anorm of the error on the finest level, see Ruede 1995.
% 
% Inputs:
% A - cell contaning system matrices
% P - cell contaning the prolongation matrices
% RHS - right hand side vector on the finest level
% numberOfLevels
% approx - approximation of A{J}^-1 RHS
% options
%   .coarseSolveName = 'backslash' compute the term coresponding to the
%       coarsest level using MATLAB backslash operator
%   .coarseSolveName = 'cg' compute the term coresponding to the coarsest
%       level using cg with a strategy for stopping described in P. Vacek
%       and J. Papež, A posteriori error estimates based on multilevel
%       decompositions (in preparation).
%       This variant requires bound on the smallest eigenvalue 
%       of the problem on the coarsest level .coarseSolveParameters.mu 
%       and parameter coarseSolveParameters.ratio in (0,1)
%   .coarseSolveName = 'none' the coarsest level term is not included in
%       the estimate
%
% Output:
% estimate of the A-norm of the error on the finest level

arguments 
    A
    P
    RHS
    J
    approx
    options.coarseSolveName = 'backslash';
    options.coarseSolveParameters = [];
end

R{J} = RHS - A{J}*approx;
z = diag(diag(A{J})) \ R{J};
estimatesq = R{J}'*z;

for j= J-1:-1:2
    R{j} = P{j+1}'*R{j+1};
    z = diag(diag(A{j})) \ R{j};   
    estimatesq = estimatesq + R{j}'*z;
end

R{1} = P{2}'*R{2};

switch(options.coarseSolveName)
    case 'none'
        estimate = sqrt(estimatesq);
    case 'backslash'
        z = options.solver.apply(A{1},R{1},zeros(size(R{1})));
        estimatesq = estimatesq + R{1}'*z;
        estimate = sqrt(estimatesq);
    case 'cg'
        [coarsestTerm,info] = cgCoarseTermMLEstimate(A{1},R{1},zeros(size(R{1})),estimatesq,options.coarseSolveParameters);
        % info.numberOfIterations
        estimatesq = estimatesq + coarsestTerm;
        estimate = sqrt(estimatesq);
end