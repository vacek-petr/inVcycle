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
    res2normRate=true,...
    numberOfVcycles=true);
results = result;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listOfRelativeTolerances = (1/2).^(1:20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with CG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for relativeTolerance = listOfRelativeTolerances

    label = ['CG (relres, \tau = ' num2str(relativeTolerance) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx);

    coarsestLevelSolver = Solver('cg');
    coarsestLevelSolver.stoppingCriterion.name = 'res2norm';
    coarsestLevelSolver.stoppingCriterion.relative = true;
    coarsestLevelSolver.stoppingCriterion.tolerance = relativeTolerance;

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approx,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        result.update(mh,approx,iter,time,coarsestLevelSolverInfo);
        result.dspl();
    end
    result.plt( ...
        marker='none', ...
        lineStyle='-',...
        errAnorm=true,...
        errAnormRate=true,...
        time=true,...
        res2norm=true,...
        res2normRate=true,...
        coarsestLevelSolverNumberOfIterations=true,...
        coarsestLevelSolverNumberOfIterationsAverage=true,...
        coarsestLevelSolverNumberOfIterationsTotal=true,...
        numberOfVcycles=true);
    results = [results result];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V-cycle with MINRES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for relativeTolerance = listOfRelativeTolerances

    label = ['MINRES (relres, \tau = ' num2str(relativeTolerance) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx);

    coarsestLevelSolver = Solver('minres');
    coarsestLevelSolver.stoppingCriterion.name = 'res2norm';
    coarsestLevelSolver.stoppingCriterion.relative = true;
    coarsestLevelSolver.stoppingCriterion.tolerance = relativeTolerance;

    while (result.errAnorm(end)/result.errAnorm(1)>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        tic
        [approx,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,coarsestLevelSolver);
        time = time + toc;

        result.update(mh,approx,iter,time,coarsestLevelSolverInfo);
        result.dspl();
    end
   result.plt( ...
        marker='none', ...
        lineStyle='-',...
        errAnorm=true,...
        errAnormRate=true,...
        time=true,...
        res2norm=true,...
        res2normRate=true,...
        coarsestLevelSolverNumberOfIterations=true,...
        coarsestLevelSolverNumberOfIterationsAverage=true,...
        coarsestLevelSolverNumberOfIterationsTotal=true,...
        numberOfVcycles=true);
    results = [results result];
end


for i=1:length(listOfRelativeTolerances)
    RESULT_CG_numberOfVcycles(i) = results(i+1).iter;
    RESULT_CG_NumberOfIterations(i) = sum(results(i+1).coarsestLevelSolverNumberOfIterations);
end
RESULT_CG_numberOfVcycles = RESULT_CG_numberOfVcycles';
RESULT_CG_NumberOfIterations = RESULT_CG_NumberOfIterations';

for j=1:length(listOfRelativeTolerances)
    RESULT_MINRES_numberOfVcycles(j) = results(j+1+length(listOfRelativeTolerances)).iter;
    RESULT_MINRES_NumberOfIterations(j) = sum(results(j+1+length(listOfRelativeTolerances)).coarsestLevelSolverNumberOfIterations);
end
RESULT_MINRES_numberOfVcycles = RESULT_MINRES_numberOfVcycles';
RESULT_MINRES_NumberOfIterations = RESULT_MINRES_NumberOfIterations';
