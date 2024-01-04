classdef Result < handle
    % RESULT class for storing the results of the analysis of V-cycle
    % method with inexact solver on the coarsest level

properties
    label
    iter % number of V-cycle iterations
    approx % matrix containing computed approximations to A{J}^(-1)*F{J}; approx(1) - initial, approx(iter+1) - last
    time % vector containing times the V-cycle iteration took
    errAnorm % vector containing A-norms of errors of the approximation after each V-cycle, errAnorm(1) - initial, errAnorm(iter+1) last
    errAnormRate % vector containing errAnorm(iter)/errAnorm(iter+1)
    err2norm % vector containing 2-norms of errors of the approximation after each V-cycle, err2norm(1) - initial, err2norm(iter+1) last
    res2norm % vector containing 2-norms of residuals of the approximation after each V-cycle, res2norm(1) - initial, res2norm(iter+1) last
    res2normRate % vector containing res2norm(iter)/res2norm(iter+1)
    
    coarsestLevelSolverNumberOfIterations % vector containing number of solver iterations in each V-cycle
    coarsestLevelSolverErrAnorm % vector containing A-norms of errors of the solver on the coarsest level in each V-cycle
    coarsestLevelSolverRes2normComputed % vector containing 2-norms of residuals of the solver on the coarsest level in each V-cycle
    coarsestLevelSolverRes2norm % vector containing 2-norms of residuals of the solver on the coarsest level (estimated inside the solver; may differ from res_2normComputed) in each V-cycle
    
    approxVcycleBackslash % matrix containing approximations computed by V-cycle with Matlab backslash on the coarsest level, approximation of x^n_ex
    gamma % vector containing (computed) parameters gamma, gamma(1) - initial, gamma(iter+1) last

    diffApproxApproxVcycleBackslashAnorm % vector containing ||x^n_ex-x^n_in||_A
    relativeDiffApproxApproxOneVcycleBackslashAnorm % vector containing ||x^new_ex - x^new_in||_A / ||x - x^prev||_A
    relativeDiffApproxApproxOneVcycleBackslash2norm % vector containing ||x^new_ex - x^new_in|| / ||x - x^prev||
    relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma % vector containing (||x^new_ex - x^new_in||_A / ||x - x^prev||_A) / gamma(end)
end
    
methods

function obj = Result(mh,label,approx,options)
    % RESULT creates result, initialize
    arguments
        mh % class MultivelFramework
        label 
        approx % initial approximation of A{J}^(-1)*F{J}
        options.gamma = []; % parameter gamma
        options.approxVcycleBackslash = []; % matrix containing approximations computed by V-cycle with Matlab backslash on the coarsest level, approximation of x^n_ex
    end
    obj.label = label;
    obj.approx = approx;
    obj.gamma = options.gamma;
    obj.approxVcycleBackslash = options.approxVcycleBackslash;

    % definitions for easier reading of the code
    J = mh.numberOfLevels;
    A = mh.A;
    F = mh.F;
    approxBackslash = mh.approxBackslash{J};
    
    obj.res2norm = norm(F{J}-A{J}*approx);
    obj.errAnorm = sqrt((approx-approxBackslash)'*A{J}*(approx-approxBackslash));
    obj.err2norm = norm(approx-approxBackslash);
end

function update(obj,mh,approx,iter,timeVcycle,coarsestLevelSolverInfo,options)
    % UPDATE update result after one V-cycle iteration
    arguments
        obj
        mh % class MatrixHierarchy
        approx % vector - computed approximation
        iter % number - current number of V-cycle iteration
        timeVcycle % number - time the Vcycle iteration took
        coarsestLevelSolverInfo % structure - information collected during the coarsest level computation
        options.approxOneVcycleBackslash = []; % vector - approximation computed using one V-cycle with Matlab backslash starting with approx(iter-1)
    end
    obj.iter = iter;
    obj.time(iter) = timeVcycle;
    obj.approx = [obj.approx approx];
    approxOneVcycleBackslash = options.approxOneVcycleBackslash;

    % definitions for easier reading of the code
    J = mh.numberOfLevels;
    A = mh.A;
    F = mh.F;
    approxBackslash = mh.approxBackslash{J};
        
    % number of solver iterations on the coarsest level
    obj.coarsestLevelSolverNumberOfIterations(iter) = coarsestLevelSolverInfo.numberOfIterations;

    % residual norms
    obj.res2norm(iter+1) = norm(F{J} - A{J}*approx);
    obj.coarsestLevelSolverRes2normComputed(iter)= coarsestLevelSolverInfo.res2normComputed;
    obj.coarsestLevelSolverRes2norm(iter) = coarsestLevelSolverInfo.res2norm;

    % error norms
    obj.coarsestLevelSolverErrAnorm(iter)= coarsestLevelSolverInfo.errAnorm;
    % mp.Digits(34);
    
    
    obj.errAnorm(iter+1) = sqrt((approx - approxBackslash)'*A{J}*(approx - approxBackslash));
    obj.err2norm(iter+1) = norm(approx - approxBackslash);

    % error and residual norms rates
    obj.errAnormRate(iter) = obj.errAnorm(iter+1)/obj.errAnorm(iter);
    obj.res2normRate(iter) = obj.res2norm(iter+1)/obj.res2norm(iter);

    % difference between approximation computed by exV-cycle and inVcycle in A norm
    if(iter <= size(obj.approxVcycleBackslash,2)-1)
        obj.diffApproxApproxVcycleBackslashAnorm(iter) = sqrt((approx - obj.approxVcycleBackslash(:,iter+1))'*A{J}*(approx - obj.approxVcycleBackslash(:,iter+1)));
    end

    % relativeDiffApproxApproxOneVcycleBackslashAnorm
    if (~isempty(approxOneVcycleBackslash))
        obj.relativeDiffApproxApproxOneVcycleBackslashAnorm(iter) = sqrt((approx-approxOneVcycleBackslash)'*A{J}*(approx-approxOneVcycleBackslash))/obj.errAnorm(iter);
        obj.relativeDiffApproxApproxOneVcycleBackslash2norm(iter) = norm(approx-approxOneVcycleBackslash)/obj.err2norm(iter);
        if(~isempty(obj.gamma))
            obj.relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma(iter) = obj.relativeDiffApproxApproxOneVcycleBackslashAnorm(iter)/obj.gamma(end);
        end
    end

end
function dspl(obj)
    % DSPL print some of the results
    messege = [ 'V-cycle iteration: ' num2str(obj.iter) newline ...
        'number of solver iterations: ' num2str(obj.coarsestLevelSolverNumberOfIterations(end)) newline ...
        'A-norm of error finest level: ' num2str(obj.errAnorm(end))  newline...
        'A-norm of error coarsest level: ' num2str(obj.coarsestLevelSolverErrAnorm(end))  newline...
        'total time: ' num2str(sum(obj.time)) newline];
    disp(messege);

end
function plt(obj,options)
    % PLT plotting the results
    arguments
        obj
        options.color = [];
        options.marker = "^";
        options.lineStyle ="-"; 
        options.errAnorm {mustBeUnderlyingType(options.errAnorm,"logical")} = false;
        options.errAnormRate {mustBeUnderlyingType(options.errAnormRate,"logical")} = false;
        options.time {mustBeUnderlyingType(options.time,"logical")} = false;
        options.res2norm {mustBeUnderlyingType(options.res2norm,"logical")} = false;
        options.res2normRate {mustBeUnderlyingType(options.res2normRate,"logical")} = false;
        options.diffApproxApproxVcycleBackslashAnorm {mustBeUnderlyingType(options.diffApproxApproxVcycleBackslashAnorm,"logical")} = false;
        options.coarsestLevelSolverNumberOfIterations {mustBeUnderlyingType(options.coarsestLevelSolverNumberOfIterations,"logical")} = false;
        options.coarsestLevelSolverNumberOfIterationsTotal {mustBeUnderlyingType(options.coarsestLevelSolverNumberOfIterationsTotal,"logical")} = false;
        options.coarsestLevelSolverNumberOfIterationsAverage {mustBeUnderlyingType(options.coarsestLevelSolverNumberOfIterationsAverage,"logical")} = false;
        options.relativeDiffApproxApproxOneVcycleBackslashAnorm {mustBeUnderlyingType(options.relativeDiffApproxApproxOneVcycleBackslashAnorm,"logical")} = false;
        options.relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma {mustBeUnderlyingType(options.relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma,"logical")} = false;
        options.relativeDiffApproxApproxOneVcycleBackslash2norm {mustBeUnderlyingType(options.relativeDiffApproxApproxOneVcycleBackslash2norm,"logical")} = false;
        options.gamma {mustBeUnderlyingType(options.gamma,"logical")} = false;
        options.numberOfVcycles {mustBeUnderlyingType(options.numberOfVcycles,"logical")} = false;
    end

    if(options.errAnorm)
        f = figure(1);
        set(f,'Name','errAnorm');
        line = semilogy(0:length(obj.errAnorm)-1, obj.errAnorm,...
            'LineWidth', 1,...
            'LineStyle',options.lineStyle,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$\| x - x^{(n)} \|_A$','interpreter','latex','FontSize',18)
    end

    if(options.errAnormRate)
        f = figure(2);
        set(f,'Name','errAnormRate');

        line = plot(1:length(obj.errAnormRate), obj.errAnormRate,...
            'LineWidth', 1,...
            'LineStyle',options.lineStyle,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$\frac{\| x - x^{(n)} \|_A}{\| x- x^{(n-1)}\|_A}$','interpreter','latex')
    end

    if (options.diffApproxApproxVcycleBackslashAnorm)
        f = figure(3);
        set(f,'Name','diffApproxApproxVcycleBackslashAnorm');
        line = semilogy(1:length(obj.diffApproxApproxVcycleBackslashAnorm), obj.diffApproxApproxVcycleBackslashAnorm,...
            'LineWidth', 1,...
            'LineStyle',options.lineStyle,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$ \| x^{(n)}_{ex} - x^{(n)}_{in} \|_A $','interpreter', 'latex')
    end

    if (options.coarsestLevelSolverNumberOfIterations)
        f = figure(4);
        set(f,'Name','numberOfSolverIterations');
        line = plot(1:length(obj.coarsestLevelSolverNumberOfIterations), obj.coarsestLevelSolverNumberOfIterations,...
            'LineWidth', 1,...
            'LineStyle',options.lineStyle,...
            'LineStyle',options.lineStyle,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('solver iter.')
    end

    if (options.coarsestLevelSolverNumberOfIterationsTotal)
        f = figure(5);
        set(f,'Name','totalNumberOfSolverIterations');
        b = bar(categorical(cellstr(obj.label)),sum(obj.coarsestLevelSolverNumberOfIterations),EdgeColor="none");
        if(options.color)
            b(1).FaceColor = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        xtips = b(1).XEndPoints;
        ytips = b(1).YEndPoints;
        labels = string(b(1).YData);
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        ylabel('total number of solver iterations');
        yticks([]);
    end

    if (options.coarsestLevelSolverNumberOfIterationsAverage)
        f = figure(6);
        set(f,'Name','averageNumberOfSolverIterations');
        b = bar(categorical(cellstr(obj.label)),mean(obj.coarsestLevelSolverNumberOfIterations),EdgeColor="none");
        if(options.color)
            b(1).FaceColor = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        xtips = b(1).XEndPoints;
        ytips = b(1).YEndPoints;
        labels = string(b(1).YData);
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        ylabel('average number of solver iterations');
    end

    if(options.time)
        f = figure(7);
        set(f,'Name','time');
        d = bar(categorical(cellstr(obj.label)),sum(obj.time),EdgeColor="none");
        if(~isempty(options.color))
            d(1).FaceColor = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        xtips = d(1).XEndPoints;
        ytips = d(1).YEndPoints;
        labels = string(d(1).YData);
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        ylabel('time');
    end

    if(options.res2norm)
        f = figure(8);
        set(f,'Name','res2norm');
        line = semilogy(0:length(obj.res2norm)-1, obj.res2norm,...
            'LineWidth', 1,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        hold on
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        legend
        xlabel('V-cycle iterations')
        ylabel('$\| b -  A x^{(n)} \|$','interpreter','latex','FontSize',18)
    end

    if(options.res2normRate)
        f = figure(9);
        set(f,'Name','res2normRate');
        line = semilogy(1:length(obj.res2normRate), obj.res2normRate,...
            'LineWidth', 1,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$\frac{\| b - A x^{(n)} \|}{\| b - Ax^{(n-1)}\|}$','interpreter','latex')
    end

    if (options.relativeDiffApproxApproxOneVcycleBackslashAnorm)
        f = figure(10);
        set(f,'Name','relativeDiffApproxApproxOneVcycleBackslashAnorm');
        line = semilogy(1:length(obj.relativeDiffApproxApproxOneVcycleBackslashAnorm), obj.relativeDiffApproxApproxOneVcycleBackslashAnorm,...
            'LineWidth', 1,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$\frac{\| x^{new}_{ex} - x^{new}_{in} \|_A}{\| x - x^{prev} \|_A}$','interpreter','latex')
    end

    if(options.relativeDiffApproxApproxOneVcycleBackslash2norm)
        f = figure(11);
        set(f,'Name','relativeDiffApproxApproxOneVcycleBackslash2norm');
        line = semilogy(1:length(obj.relativeDiffApproxApproxOneVcycleBackslash2norm), obj.relativeDiffApproxApproxOneVcycleBackslash2norm,...
            'LineWidth', 1,...
            'LineStyle',options.lineStyle,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$\frac{\| x^{new}_{ex} - x^{new}_{in} \|}{\| x - x^{prev} \|}$','interpreter','latex')
    end

    if(options.relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma)
        f = figure(12);
        set(f,'Name','relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma');
        line = semilogy(1:length(obj.relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma), obj.relativeDiffApproxApproxOneVcycleBackslashAnormScaledGamma,...
            'LineWidth', 1,...
            'LineStyle',options.lineStyle,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$\frac{\| x^{new}_{ex} - x^{new}_{in} \|_A}{\| x - x_{prev} \|_A} / \gamma $','interpreter','latex','FontSize',18)
    end

    if(options.gamma)
        f = figure(13);
        set(f,'Name','gammas');
        line = semilogy(1:length(obj.gamma), obj.gamma,...
            'LineWidth', 1,...
            'LineStyle',options.lineStyle,...
            'Marker', options.marker,...
            'DisplayName', obj.label);
        if (~isempty(options.color))
            line.Color = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        legend
        xlabel('V-cycle iterations')
        ylabel('$\gamma^{(k)}$','interpreter', 'latex')
    end

    if (options.numberOfVcycles)
        f = figure(14);
        set(f,'Name','numberOfVcycles');
        d = bar(categorical(cellstr(obj.label)),obj.iter,EdgeColor="none");
        if(~isempty(options.color))
            d(1).FaceColor = options.color;
            set(gca,'ColorOrderIndex', size(get(gca,'ColorOrder'),1) + get(gca,'ColorOrderIndex')-1);
        end
        hold on
        xtips = d(1).XEndPoints;
        ytips = d(1).YEndPoints;
        labels = string(d(1).YData);
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        yticks([]);
        ylabel('Number of V-cycles');
    end
end
end 
end
