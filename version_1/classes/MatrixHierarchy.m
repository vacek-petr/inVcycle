classdef MatrixHierarchy< handle
    % MATRIXHIERARCHY class for storing the hierarchy of stiffness and
    % prolongation matrices, their properties and the right-hand sides

    properties
        name    % string 
        A   % cell array - stiffness matrices, A{1} - coarsest level,  A{end} - finest level
        P   % cell array - prolongation matrices, P{1} empty, P{2} prolongation matrix from level 1 to level 2
        F   % cell array - right-hand sides
        numberOfLevels  % number
        numberOfUnknowns    % cell array
        ANumberOfNonZeroElements   % cell array
        ALargestEigenvalues   % cell array 
        ASmallestEigenvalues   % cell array 
        PNumberOfNonZeroElements  % cell array 
        approxBackslash  % cell array - approx of A{j}^-1F{j} computedy using Matlab backslash 
        approxBackslashTime % cell array - time of the computation with Matlab backslash
    end

    methods
        function obj = MatrixHierarchy()
        % SOLVER Construct an instance of class MatrixHierarchy
        end
        function obj = load(name)
            load([name '_A.mat'], 'A');
            load([name '_F.mat'], 'F');
            load([name '_P.mat'], 'P');
            obj.name = name;
            obj.A = A;
            obj.F = F;
            obj.P = P;
        end
        function computeProperties(obj) 
            J = size(obj.A,2);
            obj.numberOfLevels = J;
            
            A = obj.A;
            F = obj.F;
            P = obj.P;
            
            for j = 1:J
                tic
                obj.approxBackslash{j} = A{j}\F{j};
                obj.approxBackslashTime{j} = toc;
                obj.numberOfUnknowns{j} = size(A{j},2);
                obj.ANumberOfNonZeroElements{j} = nnz(A{j});
                obj.ALargestEigenvalues{j} = eigs(A{j},1,'largestabs',Tolerance=1e-4);
                obj.ASmallestEigenvalues{j} = eigs(A{j},1,'smallestabs');
            end
            for j = 2:J
                obj.PNumberOfNonZeroElements{j} = nnz(P{j});
            end
        end
        function selectLevels(obj,options)
            arguments
                obj
                options.numberOfLevels = obj.numberOfLevels;
                options.from = 1;
            end
            obj.numberOfLevels = options.numberOfLevels;
            from  = options.from;
            to = from + obj.numberOfLevels - 1;
            obj.A = obj.A(from:to);
            obj.P = obj.P(from:to);
            obj.F = obj.F(from:to);
            obj.numberOfUnknowns = obj.numberOfUnknowns(from:to);
            obj.ANumberOfNonZeroElements = obj.ANumberOfNonZeroElements(from:to);
            obj.PNumberOfNonZeroElements = obj.PNumberOfNonZeroElements(from:to);
            obj.approxBackslash = obj.approxBackslash(from:to);
            obj.approxBackslashTime = obj.approxBackslashTime(from:to);
            obj.ALargestEigenvalues = obj.ALargestEigenvalues(from:to);
            obj.ASmallestEigenvalues = obj.ASmallestEigenvalues(from:to);
        end
    end
end