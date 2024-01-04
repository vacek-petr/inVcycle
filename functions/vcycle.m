function [approx,time,coarseSolverInfo] = vcycle(A,P,j,RHS,approx,smoother,solver)
% VCYCLE Recursive implementation of the multigrid V-cycle method
% Inputs:
% A - cell array - stiffness matrices
% P - cell array - prolongation matrices
% j - number - indicating on which level in the hierarchy the function is
% RHS - vector - right hand side
% approx - vector - initial approximation of A{j}^{(-1)}RHS
% smoother - class Smoother
% solver - class Solver
% Outputs:
% approx - vector
% time - number - time the first part of the Vcycle iteration took
% coarseSolverInfo - structure containing
%   numberOfIterations
%   errAnorm
%   res2norm
%   res2normComputed
%   other outputs

if (j > 1)
    approx = smoother.apply(A{j}, RHS,approx,j,'pre');
    [correction,time,coarseSolverInfo] = vcycle(A,P,j-1,P{j}'*(RHS-A{j}*approx),zeros(size(A{j-1},1),1),smoother,solver);
    approx = smoother.apply(A{j}, RHS,approx + P{j}*correction,j,'post');
else
    [approx,coarseSolverInfo] = solver.apply(A{1}, RHS, zeros(size(RHS)));
    time = toc;

    backslashApprox = A{1}\RHS;
    coarseSolverInfo.errAnorm = sqrt((backslashApprox - approx)'*A{1}*(backslashApprox - approx));
    coarseSolverInfo.res2normComputed = norm(RHS-A{1}*approx);
    tic

end
end