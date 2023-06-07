clc; clear all; close all;

addpath('..\functions'); 
addpath('..\classes');
addpath('..\data\peak 2D BJR95');

load('peak_2D_BJR95_L8.mat','mh');
mh.selectLevels(numberOfLevels=2,from=7); % numberOfLevels=8,2, from=1,7; numberOfLelves + from has to be equal to 9
A=mh.A; P=mh.P; numberOfLevels=mh.numberOfLevels; F=mh.F{end};

vcycleRelativeTolerance = 10^(-11);
vcycleMaximumNumberOfIterations = 100;
smoother = Smoother('gs',0,3);

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

while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
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

listOfGammas = (1/2).^(1:20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG GRML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['CG (GRML, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'GRML';
    coarsestLevelSolver.stoppingCriterion.mu = mh.ASmallestEigenvalues{1};
    coarsestLevelSolver.stoppingCriterion.gamma = gamma;
    coarsestLevelSolver.stoppingCriterion.cs = 0.85;
    coarsestLevelSolver.stoppingCriterion.fineLevelsContribution = computeErrAnormMLEstimate(mh.A,mh.P,mh.F{end},mh.numberOfLevels,approx,coarseSolveName='none')^2;

    coarseSolveParameters.ratio = 0.1;
    coarseSolveParameters.mu = mh.ASmallestEigenvalues{1};
    errAnormMLEstimate = computeErrAnormMLEstimate(mh.A,mh.P,mh.F{end},mh.numberOfLevels,approx,coarseSolveName='cg',coarseSolveParameters=coarseSolveParameters);
    
    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        coarsestLevelSolver.stoppingCriterion.fineLevelsContribution = computeErrAnormMLEstimate(mh.A,mh.P,mh.F{end},mh.numberOfLevels,approx,coarseSolveName='none')^2;
        errAnormMLEstimateNew = computeErrAnormMLEstimate(mh.A,mh.P,mh.F{end},mh.numberOfLevels,approx,coarseSolveName='cg',coarseSolveParameters=coarseSolveParameters);
        coarsestLevelSolver.stoppingCriterion.gamma = (2^(1/(10))-1)*(coarsestLevelSolver.stoppingCriterion.cs)*(errAnormMLEstimateNew/errAnormMLEstimate);
        result.gamma(iter+1) =coarsestLevelSolver.stoppingCriterion.gamma;

        errAnormMLEstimate = errAnormMLEstimateNew;
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
        relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma=true, ...
        gamma=true);
    results = [results result];
end

for i=1:length(listOfGammas)
    RESULT_CG_numberOfVcycles(i) = results(i+1).iter;
    RESULT_CG_NumberOfIterations(i) = sum(results(i+1).coarsestLevelSolverNumberOfIterations);
end
RESULT_CG_numberOfVcycles = RESULT_CG_numberOfVcycles';
RESULT_CG_NumberOfIterations = RESULT_CG_NumberOfIterations';