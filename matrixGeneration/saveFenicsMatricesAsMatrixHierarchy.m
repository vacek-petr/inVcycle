% script for saving matrices generated in Fenics as MatrixHierarchy class
% object

% load "A.mat", "F.mat", "P.mat", and "BN.mat"

addpath('..\classes');

mh = MatrixHierarchy();
mh.name = 'Poisson'; % or jump-1024 
mh.A = A;
J = size(mh.A,2);
mh.F = F;
mh.P = {};
mh.P(2:J) = P;

% cut rows and collums coresponding to nodes on boundary
% this is done in order to have matrices satisfying the Galerkin condition
for j=1:J
    mh.A{j}(BN{j},:)=[];
    mh.A{j}(:,BN{j})=[];
    mh.F{j}(BN{j})=[];
    mh.F{j} = mh.F{j}';
    if(j ~= 1)
        mh.P{j}(:,BN{j-1})=[];
        mh.P{j}(BN{j},:)=[];
        mh.P{j} = mh.P{j}.*(mh.P{j} ~= 0); % drop zero entries from sparse matrix
    end
end
mh.computeProperties();
save(mh.name,"mh")