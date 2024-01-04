clc; clear all; close all;

addpath('..\functions'); 
addpath('..\classes');
addpath('..\data\Poisson\')
addpath('..\data\jump-1024\')
load('Poisson.mat','mh'); % or load('jump-1024.mat','mh');

mh.selectLevels(numberOfLevels=6,from=2);

% mh.switchToMp(34); turn on simulated quad precision using Advapix toolbox

A=mh.A; P=mh.P; numberOfLevels=mh.numberOfLevels; F=mh.F{end};

vcycleAbsoluteTolerance = 10^(-11);
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

listOfGammas = [1e-4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate T2norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T2norm = estimateT2norm(mh); 
% T2norm = 32.5157; % jump-1024
T2norm = 31.8365; % Poisson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG relres %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['CG (relres, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'res2norm';
    coarsestLevelSolver.stoppingCriterion.relative = true;
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(mh.ASmallestEigenvalues{1})*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*1/T2norm;

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
    result.errAnorm = double(result.errAnorm');
    result.errAnormRate = double(result.errAnormRate');
    result.relativeDiffApproxApproxOneVcycleBackslashAnorm= double(result.relativeDiffApproxApproxOneVcycleBackslashAnorm');
    result.coarsestLevelSolverNumberOfIterations = double(result.coarsestLevelSolverNumberOfIterations');
end
