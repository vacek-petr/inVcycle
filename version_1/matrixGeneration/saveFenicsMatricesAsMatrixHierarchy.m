% script for saving matrices generated in Fenics as MatrixHierarchy class
% object

addpath('data_new\peak 2D BJR95'); addpath('classes');

mh = MatrixHierarchy();
mh.name = 'peak_2D_BJR95_L8';
load("A.mat");
mh.A = A;
J = size(mh.A,2);
load("F.mat");
mh.F = F;
load("P.mat");
mh.P = {};
mh.P(2:J) = P;
load("BN.mat")

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