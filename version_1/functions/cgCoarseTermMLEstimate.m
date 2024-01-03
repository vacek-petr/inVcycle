function [sumAlpha2normRes,info] = cgCoarseTermMLEstimate(A,b,approx,fineLevelsContribution,parameters)
% CGCOARSETERMMLESTIMATE Implementation of the CG algorithm for computing the approximation
% of b'*A^(-1)*b, in the context of multilevel estimator (Ruede 1995).
% The computation is stopped using the estimate on the A-norm of the error based on the Gauss-Radau
% quadrature (Meurant,Tichy 2019, equation(6)), when 
% ||error||^2_A < (fineLevelsContribution + sumAlpha2normRes)*ratio
% see J. Papež and P. Vacek A posteriori error estimates based on 
% multilevel decompositions(in preparation).
% 
% Inputs:
% A - symmetric positive definite matrix
% b - right-hand side vector
% approx - initial approximation of A^(-1)*b
% parameters - structure - must contain
%   mu - lower estimate of the smallest eigenvalue of A
%   ratio
% fineLevelsContribution of the MLestimator
%   
% Outputs:
% sumAlpha2normRes- approximation of b'*A^(-1)*b (the term in ML estimator corresponding to the coarsest level)
% info - structure containing
%   res2norm  - Euclidean norm of the residual b-A*approx
%   numberOfIterations

iter = 0;
r = b - A*approx;
rr_new = r'*r;
p = r;

mu = parameters.mu;
gamma_mu = 1/mu;
ratio = parameters.ratio;
sumAlpha2normRes = 0;
UpperBoundGRsq = gamma_mu*rr_new;
if(UpperBoundGRsq > ratio*fineLevelsContribution)
    stoppingCriterionNotSatisfied = 1;
else
    stoppingCriterionNotSatisfied = 0;
end

while(stoppingCriterionNotSatisfied)

    rr = rr_new;
    iter = iter + 1;
    Ap = A*p;
    alpha = rr/(p'*Ap);
    approx = approx + alpha*p;
    r = r - alpha*Ap;
    rr_new = r'*r;
    delta = rr_new/rr;
    p = r + delta*p;

    gamma_mu = (gamma_mu-alpha)/(mu*(gamma_mu-alpha) + delta);
    UpperBoundGRsq = gamma_mu*rr_new;
    sumAlpha2normRes = sumAlpha2normRes + alpha*rr;
    if (UpperBoundGRsq > ratio*(fineLevelsContribution + sumAlpha2normRes))
        stoppingCriterionNotSatisfied = 1;
    else
        stoppingCriterionNotSatisfied = 0;
    end
end
info.res2norm = sqrt(rr_new);
info.numberOfIterations = iter;
end
