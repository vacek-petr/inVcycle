classdef Solver
% SOLVER class for calling solver for Ax = b

properties
    name % string
    stoppingCriterion % structure
    parameters % structure
end

methods
    function obj = Solver(name,options)
    % SOLVER Construct an instance of class Solver
    arguments
        name {mustBeMember(name,{'cg','minres','backslash','blr'})}
        options.stoppingCriterion = []
        options.parameters = []
    end
    obj.name = name;
    obj.stoppingCriterion = options.stoppingCriterion;
    obj.parameters = options.parameters;
end

function [approx,info] = apply(obj,A,b,approx)
    switch(obj.name)
        case 'cg'
            [approx,info] = cg(A,b,approx,obj.stoppingCriterion);
        case 'minres'
            [approx,info] = my_minres(A,b,approx,obj.stoppingCriterion);
        case 'backslash'
            [approx,info] = backslash(A,b);
        case 'blr'
            [approx,info] = blr(A,b,obj.parameters.BLRLU);
    end
end
end
end