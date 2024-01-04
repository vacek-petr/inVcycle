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
%
% Outputs:
% approx - approximation of A^(-1)*b
% info - structure containing
%   res2norm - Euclidean norm of the residual b-A*approx
%   numberOfIterations - number of performed CG iterations

% delete

plotErrAnormAndBounds = 1;

if (plotErrAnormAndBounds)
    errAnormPlot = [];
    upperBoundGRPlot = [];
end

iter = 0;
r = b - A*approx;
rr_new = r'*r;
p = r;

stoppingCriterionName = stoppingCriterion.name;

switch(stoppingCriterionName)
    case 'errAnorm'
        tolerance = stoppingCriterion.tolerance;
        approxBackslash = A\b;
        errAnorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
        if (errAnorm > tolerance)
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
        if (plotErrAnormAndBounds)
            approxBackslash = A\b; 
        end
        tolerance = stoppingCriterion.tolerance;
        mu = stoppingCriterion.mu;
        gamma_mu = 1/mu;
        UpperBoundGRsq = gamma_mu*rr_new;
        if (UpperBoundGRsq < 0)
            disp('GR bound failed!');
            pause;
        end
        if(sqrt(UpperBoundGRsq) > tolerance)
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
            errAnorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
            if (plotErrAnormAndBounds)
                errAnormPlot(iter) = errAnorm;
            end
            if (errAnorm > tolerance)
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
            if (plotErrAnormAndBounds)
                errAnorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
                errAnormPlot(iter) = errAnorm;
            end
            gamma_mu = (gamma_mu-alpha)/(mu*(gamma_mu-alpha) + delta);
            UpperBoundGRsq = gamma_mu*rr_new;
            upperBoundGRPlot(iter) = UpperBoundGRsq;
            if (UpperBoundGRsq < 0)
                disp('GR bound failed.');
                pause;
            end
            if (sqrt(UpperBoundGRsq) > tolerance)
                stoppingCriterionNotSatisfied = 1;
            else
                stoppingCriterionNotSatisfied = 0;
            end
    end
end
if (plotErrAnormAndBounds)
    switch(stoppingCriterionName)
        case 'errAnorm'
                f = figure(15);
                set(f,'Name','CG errAnorm and bound');
                semilogy(1:iter,errAnormPlot)
                hold on
                plot(1:iter,ones(iter,1)*tolerance)
                hold on
        case 'GR'
                f = figure(15);
                set(f,'Name','CG errAnorm and bound');
                semilogy(1:iter,errAnormPlot)
                hold on
                plot(1:iter,ones(iter,1)*tolerance)
                hold on
                semilogy(1:iter,sqrt(upperBoundGRPlot))
                hold on
    end
end

info.res2norm = sqrt(rr_new);
info.numberOfIterations = iter;
end