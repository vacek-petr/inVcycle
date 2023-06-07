function T2norm = estimateT2norm(mh)
% T2NORM function for estimating the Euclidean norm of the restriction
% operator T = P{2}'*P{3}'*...*P{end}'

P = mh.P{2};
T = P';

for j=3:mh.numberOfLevels
    P = mh.P{j};
    T = T*P';
end
T2norm = normest(T,1e-4);

end
