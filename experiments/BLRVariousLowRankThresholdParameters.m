% The BLR LU decompositions of system matrices are not part of the repository.
% They can be generate by running the script computeBLRLU.m

clc; clear all; close all;

addpath('..\functions'); 
addpath('..\classes');
addpath('..\data\peak 2D BJR95'); 
addpath('..\data\peak 2D BJR95\BLR LU decompositions\');

load('peak_2D_BJR95_L8.mat','mh');
numberOfLevels = 6; % 6,5,4; numberOfLelves + from has to be equal to 9
from = 3; % 3,4,5 
mh.selectLevels(numberOfLevels=numberOfLevels,from=from);
A=mh.A; P=mh.P; F=mh.F{end};

vcycleRelativeTolerance = 10^(-11);
vcycleMaximumNumberOfIterations = 30;
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
    iter = iter+1;
    
    tic
    [approx,time,coarsestLevelSolverInfo] = vcycle(A,P,numberOfLevels,F,approx,smoother,Solver("backslash"));                         
    time = time+toc;
    
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
% V-cycle with BLR direct solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load BLR LU decomposition computed by function computeBLRLU.m
BLRLUFileNameFirstPart = ['LU_J' num2str(from) '_uprec_double_lrth'];

for i = 1:24
    BLRLUFileName = [BLRLUFileNameFirstPart num2str(i) '.mat'];
    load(BLRLUFileName);

    label = ['BLR (\epsilon='  num2str(0.5^i) ')'];
    disp(['Using V-cycle with ' label  newline]);

    iter = 0;
    approx = initialApprox;
    result = Result(mh,label,approx,approxVcycleBackslash = results(1).approx);

    coarsestLevelSolver = Solver('blr');
    coarsestLevelSolver.parameters.BLRLU = BLRLU;
    
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
        relativeDiffApproxApproxOneVcycleBackslashAnorm=true,...
        numberOfVcycles=true);
    results = [results result];
end

for j = 1:length(results)-1
    RESULT_BLR_numberOfVcycles(j) = results(j+1).iter;
end
RESULT_BLR_numberOfVcycles = RESULT_BLR_numberOfVcycles';
