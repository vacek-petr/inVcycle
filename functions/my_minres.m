function [approx,info] = my_minres(A,b,approx,stoppingCriterion)
% MY_MINRES method derived in (Paige, Saunders 1975). This
% implementation is based on the algorithm in 
% (Sleijpen, Van der Vorst, Modersitzki 2001). We consider various 
% stopping various stopping criteria.
%
% Inputs:
% A - matrix
% b - right-hand side vector
% approx - initial approximation of A^(-1)*b
% stoppingCriterion - structure
%   implemented variants:
%   .name='errAnorm' - stopping when (approximate) absolute A norm of the error is less than given .tolerance
%   .name='res2norm' - stopping when (relative [.relative=1] or absolute [.relative=0]) Euclidien norm of residual is less
%       than given .tolerance 
%
% Outputs:
% approx - approximation of A^(-1)*b
% info - structure containing
%   res2norm - Euclidean norm of the residual b-A*approx
%   numberOfIterations - number of performed CG iterations

iter = 0;
r = b - A*approx;
r2norm = norm(r);
v = r/r2norm;
beta = 0;
betaTilde = 0;
c = -1;
s = 0;
vOld = zeros(size(approx));
w = zeros(size(approx));
wDobleTilde = v;

switch(stoppingCriterion.name)
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
        if(r2norm > tolerance)
            stoppingCriterionNotSatisfied = 1;
        else
            stoppingCriterionNotSatisfied = 0;
        end
end


while(stoppingCriterionNotSatisfied)
    
    iter = iter + 1;
    vTilde = A*v - beta*vOld;
    alpha = v'*vTilde;
    vTilde = vTilde - alpha*v;
    beta = norm(vTilde);
    vOld = v;
    v = vTilde/beta;
    l1 = s*alpha - c*betaTilde;
    l2 = s*beta;
    alphaTilde = -s*betaTilde-c*alpha;
    betaTilde = c*beta;
    l0 = sqrt(alphaTilde^2 + beta^2);
    c = alphaTilde/l0;
    s = beta/l0;
    wTilde = wDobleTilde - l1*w;
    wDobleTilde = v - l2*w;
    w = wTilde/l0;
    approx = approx + r2norm*c*w;
    r2norm = s*r2norm;

    switch(stoppingCriterion.name)
        case 'errAnorm'
            err_Anorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
            if (err_Anorm > tolerance)
                stoppingCriterionNotSatisfied = 1;
            else
                stoppingCriterionNotSatisfied = 0;
            end
        case 'res2norm'
            if(r2norm > tolerance)
                stoppingCriterionNotSatisfied = 1;
            else
                stoppingCriterionNotSatisfied = 0;
            end
    end
end
    info.numberOfIterations = iter;
    info.res2norm = r2norm;
end