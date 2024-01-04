clc; clear all; close all;

addpath('..\functions'); 
addpath('..\classes');
addpath('..\data\Poisson\')
addpath('..\data\jump-1024\')
load('Poisson.mat','mh'); % or load('jump-1024.mat','mh');

mh.selectLevels(numberOfLevels=6,from=2);

A=mh.A; P=mh.P; numberOfLevels=mh.numberOfLevels; F=mh.F{end};

vcycleAbsoluteTolerance = 0; %  10^(-11) or 0 (if we want to stop after vcycleMaximumNumberOfIterations)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listOfGammas = [0.3,1e-3,1e-4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG errERR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['CG (errERR, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash=results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'errAnorm';
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*result.errAnorm(end);

    while (result.errAnorm(end)>vcycleAbsoluteTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;

        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);
        result.dspl();

        coarsestLevelSolver.stoppingCriterion.tolerance = gamma*result.errAnorm(end);
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
    result.errAnormRate = result.errAnormRate';
    result.relativeDiffApproxApproxOneVcycleBackslashAnorm = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';
    result.coarsestLevelSolverNumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
end
