function [approx,info] = cg(A,b,approx,stoppingCriterion)
% CG Implementation of the CG algorithm (Hestenes, Stiefel 1952) with
% various stopping criteria.
%
% Inputs:
% A - symmetric positive definite matrix
% b - right-hand side vector
% approx - initial approximation of A^(-1)*b
% stoppingCriterion - structure
%   implemented variants:
%   .name='errAnorm' - stopping when (approximate) absolute A norm of the error is less than given .tolerance
%   .name='res2norm' - stopping when (relative [.relative=1] or absolute [.relative=0]) Euclidien norm of residual is less
%       than given .tolerance
%   .name='GR' - stopping criterion based on Gauss-Radau quadrature
%       interpretation of CG, see (Meurant,Tichy 2019, equation(6)),
%       requires .tolerance and .mu - lower bound on the smallest eigenvalue of A 
%   .name='GRML' - stopping criterion for CG when used as coarsest-level
%       solver in multigrid V-cycle, see (Vacek, Carson, Soodhalter, The
%       effect of approximate coarsest-level solves...) requires
%       .cs - constant in the ML estimator lower bound, 
%       .fineLevelsContribution - of the ML estimator (squared)
%       .gamma 
%
% Outputs:
% approx - approximation of A^(-1)*b
% info - structure containing
%   res2norm - Euclidean norm of the residual b-A*approx
%   numberOfIterations - number of performed CG iterations


iter = 0;
r = b - A*approx;
rr_new = r'*r;
p = r;

stoppingCriterionName = stoppingCriterion.name;

switch(stoppingCriterionName)
    case 'errAnorm'
        tolerance = stoppingCriterion.tolerance;
        approxBackslash = A\b;
        err_Anorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
        if (err_Anorm > tolerance)
            stoppingCriterionNotSatisfied = 1;
        else
            stoppingCriterionNotSatisfied = 0;
        end
    case 'res2norm'
        if (stoppingCriterion.relative)
            tolerance = norm(b)*stoppingCriterion.tolerance;
        else
            tolerance = stoppingCriterion.tolerance;
        end
        if(sqrt(rr_new) > tolerance)
            stoppingCriterionNotSatisfied = 1;
        else
            stoppingCriterionNotSatisfied = 0;
        end
    case 'GR'
        tolerance = stoppingCriterion.tolerance;
        mu = stoppingCriterion.mu;
        gamma_mu = 1/mu;
        UpperBoundGRsq = gamma_mu*rr_new;
        if(sqrt(UpperBoundGRsq) > tolerance)
            stoppingCriterionNotSatisfied = 1;
        else
            stoppingCriterionNotSatisfied = 0;
        end    
    case 'GRML'
        fineLevelsContribution = stoppingCriterion.fineLevelsContribution;
        gamma = stoppingCriterion.gamma;
        cs = stoppingCriterion.cs;
        mu = stoppingCriterion.mu;
        gamma_mu = 1/mu;
        sumAlpha2normRes = 0;
        UpperBoundGRsq = gamma_mu*rr_new;
        if(UpperBoundGRsq > cs*gamma^2/(1-cs*gamma^2)*(fineLevelsContribution))
            stoppingCriterionNotSatisfied = 1;
        else
            stoppingCriterionNotSatisfied = 0;
        end
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

    switch(stoppingCriterionName)
        case 'errAnorm'
            err_Anorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
            if (err_Anorm > tolerance)
                stoppingCriterionNotSatisfied = 1;
            else
                stoppingCriterionNotSatisfied = 0;
            end
        case 'res2norm'
            if(sqrt(rr_new) > tolerance)
                stoppingCriterionNotSatisfied = 1;
            else
                stoppingCriterionNotSatisfied = 0;
            end
        case 'GR'
            gamma_mu = (gamma_mu-alpha)/(mu*(gamma_mu-alpha) + delta);
            UpperBoundGRsq = gamma_mu*rr_new;
            if (sqrt(UpperBoundGRsq) > tolerance)
                stoppingCriterionNotSatisfied = 1;
            else
                stoppingCriterionNotSatisfied = 0;
            end
        case 'GRML'
            gamma_mu = (gamma_mu-alpha)/(mu*(gamma_mu-alpha) + delta);
            UpperBoundGRsq = gamma_mu*rr_new;
            sumAlpha2normRes = sumAlpha2normRes + alpha*rr;
            if (UpperBoundGRsq > cs*gamma^2/(1-cs*gamma^2)*(fineLevelsContribution + sumAlpha2normRes))
                stoppingCriterionNotSatisfied = 1;
            else
                stoppingCriterionNotSatisfied = 0;
            end
    end
end
info.res2norm = sqrt(rr_new);
info.numberOfIterations = iter;
end