clc; clear all; close all;

addpath('..\functions'); 
addpath('..\classes');
addpath('..\data\peak 2D BJR95');

load('peak_2D_BJR95_L8.mat','mh');
mh.selectLevels(numberOfLevels=8,from=1); % numberOfLevels=8,2, from=1,7; numberOfLelves + from has to be equal to 9
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

listOfGammas = (1/2).^(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG errerrr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['CG (errerr, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash=results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'errAnorm';
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*result.errAnorm(end);

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        coarsestLevelSolver.stoppingCriterion.tolerance = gamma*result.errAnorm(end);
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

RESULT_CG_errerr_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_CG_errerr_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with MINRES errerr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['MINRES (errerr, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('minres');
    coarsestLevelSolver.stoppingCriterion.name = 'errAnorm';
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*result.errAnorm(end);

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        coarsestLevelSolver.stoppingCriterion.tolerance = gamma*result.errAnorm(end);
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

RESULT_MINRES_errerr_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_MINRES_errerr_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate T2norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T2norm = estimateT2norm(mh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
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

RESULT_CG_relres_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_CG_relres_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with MINRES relres %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['MINRES (relres, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('minres');
    coarsestLevelSolver.stoppingCriterion.name = 'res2norm';
    coarsestLevelSolver.stoppingCriterion.relative = true;
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(mh.ASmallestEigenvalues{1})*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*1/T2norm;

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
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

RESULT_MINRES_relres_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_MINRES_relres_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG resres %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['CG (resres, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'res2norm';
    coarsestLevelSolver.stoppingCriterion.relative = false;
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(mh.ASmallestEigenvalues{1})*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*result.res2norm(end);

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

       approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(mh.ASmallestEigenvalues{1})*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*result.res2norm(end);
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

RESULT_CG_resres_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_CG_resres_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with MINRES resres %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['MINRES (resres, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('minres');
    coarsestLevelSolver.stoppingCriterion.name = 'res2norm';
    coarsestLevelSolver.stoppingCriterion.relative = false;
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(mh.ASmallestEigenvalues{1})*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*result.res2norm(end);

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

       approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(mh.ASmallestEigenvalues{1})*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*result.res2norm(end);
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

RESULT_MINRES_resres_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_MINRES_resres_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG GRML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        coarsestLevelSolver.stoppingCriterion.fineLevelsContribution = computeErrAnormMLEstimate(mh.A,mh.P,mh.F{end},mh.numberOfLevels,approx,coarseSolveName='none')^2;
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

RESULT_CG_GRML_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_CG_GRML_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG GRres %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for gamma = listOfGammas

    label = ['CG (GRres, \gamma = ' num2str(gamma) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx,gamma=gamma);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'GR';
    coarsestLevelSolver.stoppingCriterion.mu = mh.ASmallestEigenvalues{1};
    coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*result.res2norm(end);

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approxOneVcycle,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

       approxOneVcycleBackslash = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));
        approx = approxOneVcycle;
        result.update(mh,approx,iter,time,coarsestLevelSolverInfo,approxOneVcycleBackslash=approxOneVcycleBackslash);

        coarsestLevelSolver.stoppingCriterion.tolerance = gamma*sqrt(1/mh.ALargestEigenvalues{mh.numberOfLevels})*result.res2norm(end);
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

RESULT_CG_GRres_NumberOfIterations = result.coarsestLevelSolverNumberOfIterations';
RESULT_CG_GRres_relativeDiff = result.relativeDiffApproxApproxOneVcycleBackslashAnorm';
