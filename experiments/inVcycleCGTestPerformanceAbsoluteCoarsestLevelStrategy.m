clc; clear all; close all;

addpath('..\functions'); 
addpath('..\classes');
addpath('..\data\Poisson\')
addpath('..\data\jump-1024\')
load('Poisson.mat','mh'); % or load('Poisson.mat','mh');

mh.selectLevels(numberOfLevels=6,from=2);  % numberOfLevels= 6 or 3 ,from= 2 or 5

A=mh.A; P=mh.P; numberOfLevels=mh.numberOfLevels; F=mh.F{end};

vcycleAbsoluteTolerance = 10^(-11); % or 10^(-4);
vcycleMaximumNumberOfIterations = 50;
smoother = Smoother('gs_sym',1,1);

initialApprox = zeros(size(F));

disp('Using Matlab backslash on finest level');
disp(['solution time: ' num2str(mh.approxBackslashTime{end}) newline]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with Matlab backslash %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
label = 'Matlab backslash';
disp(['Using V-cycle with ' label  newline]);

iter = 0;
approx = initialApprox;
result = Result(mh,label,approx);

while (result.errAnorm(end)>vcycleAbsoluteTolerance)&&(iter<vcycleMaximumNumberOfIterations)
    iter = iter + 1;
    
    tic
    [approx,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));                         
    time = time + toc;
    
    result.update(mh,approx,iter,time,coarsestLevelSolverInfo);
    result.dspl();
end

result.plt( ...
    marker='none', ...
    color='black',...
    lineStyle='-',...
    errAnorm=true,...
    errAnormRate=true,...
    time=true,...
    res2norm=true,...
    res2normRate=true);
results = result;

RESULT_BACKSLASH_errAnorm = result.errAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set epsilon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listOfEpsilons = (1-2/3)*vcycleAbsoluteTolerance; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG abs err %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for epsilon = listOfEpsilons

    label = ['CG (errAnorm, \epsilon = ' num2str(epsilon) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash=results(1).approx,gamma=epsilon);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'errAnorm';
    coarsestLevelSolver.stoppingCriterion.tolerance = epsilon;

    while (result.errAnorm(end)>vcycleAbsoluteTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);        
        result.dspl();
    end
    result.plt( ...
        marker='none',...
        lineStyle='-',...
        errAnorm=true,...
        errAnormRate=true,...
        time=true,...
        res2norm=true,...
        res2normRate=true,...
        coarsestLevelSolverNumberOfIterations=true,...
        coarsestLevelSolverNumberOfIterationsAverage=true,...
        coarsestLevelSolverNumberOfIterationsTotal=true,...
        numberOfVcycles=true, ...
        relativeDiffApproxApproxOneVcycleBackslash2norm=true, ...
        diffApproxApproxVcycleBackslashAnorm=true, ...
        relativeDiffApproxApproxOneVcycleBackslashAnorm=true,...
        relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma=true);
    results = [results result];
end

RESULT_CG_ERR_errAnorm = result.errAnorm';
RESULT_CG_ERR_numberOfIter = result.coarsestLevelSolverNumberOfIterations';
RESULT_CG_ERR_diff = result.diffApproxApproxVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG abs GR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for epsilon = listOfEpsilons

    label = ['CG (GR, \epsilon = ' num2str(epsilon) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash=results(1).approx,gamma=epsilon);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'GR';
    coarsestLevelSolver.stoppingCriterion.mu = (1-10^(-3))*mh.ASmallestEigenvalues{1}; % factor (1-10^(-3)) insure that mu is a lower bar on the smallest eigenvalue
    coarsestLevelSolver.stoppingCriterion.tolerance = epsilon;

    while (result.errAnorm(end)>vcycleAbsoluteTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        result.dspl();
    end
    result.plt( ...
        marker='none',...
        lineStyle='-',...
        errAnorm=true,...
        errAnormRate=true,...
        time=true,...
        res2norm=true,...
        res2normRate=true,...
        coarsestLevelSolverNumberOfIterations=true,...
        coarsestLevelSolverNumberOfIterationsAverage=true,...
        coarsestLevelSolverNumberOfIterationsTotal=true,...
        numberOfVcycles=true, ...
        relativeDiffApproxApproxOneVcycleBackslash2norm=true, ...
        diffApproxApproxVcycleBackslashAnorm=true, ...
        relativeDiffApproxApproxOneVcycleBackslashAnorm=true,...
        relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma=true);
    results = [results result];
end

RESULT_CG_GR_errAnorm = result.errAnorm';
RESULT_CG_GR_diff = result.diffApproxApproxVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG res %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for epsilon = listOfEpsilons

    label = ['CG (res, \epsilon = ' num2str(epsilon) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash=results(1).approx,gamma=epsilon);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'res2norm';
    coarsestLevelSolver.stoppingCriterion.relative = false;
    coarsestLevelSolver.stoppingCriterion.tolerance = epsilon*sqrt(mh.ASmallestEigenvalues{1});

    while (result.errAnorm(end)>vcycleAbsoluteTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);
        result.dspl();
    end
    result.plt( ...
        marker='none',...
        lineStyle='-',...
        errAnorm=true,...
        errAnormRate=true,...
        time=true,...
        res2norm=true,...
        res2normRate=true,...
        coarsestLevelSolverNumberOfIterations=true,...
        coarsestLevelSolverNumberOfIterationsAverage=true,...
        coarsestLevelSolverNumberOfIterationsTotal=true,...
        numberOfVcycles=true, ...
        relativeDiffApproxApproxOneVcycleBackslash2norm=true, ...
        diffApproxApproxVcycleBackslashAnorm=true, ...
        relativeDiffApproxApproxOneVcycleBackslashAnorm=true,...
        relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma=true);
    results = [results result];
end


for result = results
    result.errAnorm = result.errAnorm';
    result.res2norm = result.res2norm';
    result.diffApproxApproxVcycleBackslashAnorm= result.diffApproxApproxVcycleBackslashAnorm';
    result.coarsestLevelSolverNumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
end