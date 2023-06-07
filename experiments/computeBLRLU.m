% script for computing the BLR LU decompositions of coarsest-level matrices
% used in the blr solver
%
% 
% modification of the code at
% https://gitlab.com/theo.andreas.mary/BLRstability
% (Higham, Mary 2021)
%
% Possible parameter settings:
% betatype: threshold type (beta parameter)
%     1 = local
%     2 = local scaled
%     3 = global
%     4 = global scaled
%   precision: floating-point precision
%     'd' = double  - most of the computation done in sparse, to compute the
%           SVD of the block we switch to dense and right afterwards back
%           to sparse
%     's' = single  - sparse matrices cast to dense, resulting BLR LU is cast
%           back to sparse double
%     'h' = half  - sparse matrices cast to dense single, SVD done in single,
%           matmul implemented using chop (half), resulting BLR LU is cast
%           back to sparse double. This setting requires chop.m from 
%           the Matrix Computations Toolbox 
%           http://www.ma.man.ac.uk/~higham/mctoolbox/
%   variant: BLR factorization algorithm
%     'UFC' = Update-Factor-Compress
%     'UCF' = Update-Compress-Factor
%     'FUC' = Factor-Update-Compress
%     'CUF' = Compress-Update-Factor
%   recomp: intermediate recompressions (0 = off, 1 = on)
%   terminalOutputFileName: output file where my_fprintf redirects  (if '', no redirection)

clear all; close all; clc;

load('peak_2D_BJR95_L8.mat','mh');
Afull = mh.A{5}; % 3,4,5
parameters.blockSize = 53; % 39,79,53
parameters.betatype = 4;
parameters.variant = 'UCF';
parameters.recomp = 1;
parameters.precision = 'h'; % 'd','s','h'

BLRLUFileName = 'LU_J5_uprec_half_lrth';
terminalOutputFileName = '';

if strcmp(terminalOutputFileName,'')
    fout = 0;
else
    fout = fopen(terminalOutputFileName,'a');
end

switch(parameters.precision)
    case 'd', precisionStr = 'double';
    case 's', precisionStr = 'single';
    case 'h', precisionStr = 'half';
    otherwise, error('Wrong uprec parameter\n');
end

switch(parameters.betatype)
    case 1, betatypeStr = 'local';
    case 2, betatypeStr = 'local, scaled';
    case 3, betatypeStr = 'global';
    case 4, betatypeStr = 'global, scaled';
    otherwise, error('Wrong betatype parameter\n');
end

A = constructBlockMatrix(Afull,parameters.blockSize);

parameters.numberOfBlocks = size(A,1);
if parameters.precision=='d'
    Au = A;
else
    Au = cell(parameters.numberOfBlocks);
    for i=1:parameters.numberOfBlocks
        for j=1:parameters.numberOfBlocks
            Au{i,j} = single(full(A{i,j}));
        end
    end
end


for i = 1:24
    parameters.lowRankThreshold = 1/2^(i);

    my_fprintf(fout,{'====================================\n'});
    my_fprintf(fout,{'Matrix of order           n = %d\n',size(Afull,1)});
    my_fprintf(fout,{'Number of block-rows/cols p = %d\n',parameters.numberOfBlocks});
    my_fprintf(fout,{'Block size                b = %d\n',parameters.blockSize });
    my_fprintf(fout,{'Low-rank threshold          = %.0e\n',parameters.lowRankThreshold});
    my_fprintf(fout,{'Threshold type (beta param) = %s\n',betatypeStr});
    my_fprintf(fout,{'Floating-point precision    = %s\n',precisionStr});
    my_fprintf(fout,{'BLR factorization algorithm = %s\n',parameters.variant});
    my_fprintf(fout,{'Intermediate recompressions = %d\n',parameters.recomp});
    my_fprintf(fout,{'====================================\n'});

    BLRLU = blrlu(Au,Afull,parameters);

    if parameters.precision ~='d'
        for k=1:parameters.numberOfBlocks
            for j=1:parameters.numberOfBlocks
                BLRLU{k,j}.X = sparse(double(BLRLU{k,j}.X)); % switch back to sparse double
                if BLRLU{k,j}.IsLR
                    BLRLU{k,j}.Y = sparse(double(BLRLU{k,j}.Y)); % switch back to sparse double
                end
            end
        end
    end


    err = bwderr(Afull,BLRLU,parameters.numberOfBlocks);
    my_fprintf(fout,{'BLR LU backward error       = %.2e\n',err});
    my_fprintf(fout,{'====================================\n'});

    filename = [BLRLUFileName num2str(i) '.mat' ];
    save(filename,'BLRLU')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = constructBlockMatrix(Afull,blksiz)
n = size(Afull,1);
clusters = 1:blksiz:n+1; % vector contaning the starts of the clusters
blknum = length(clusters)-1;
A = cell(blknum);
for i=1:blknum
    for j=1:blknum
        A{i,j} = Afull(clusters(i):clusters(i+1)-1,clusters(j):clusters(j+1)-1);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n = mynorm(A)
n = norm(A,'fro');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LU = blrlu(A,Afull,parameters)
% computes BLR LU factors LU = A
% LU.X = row basis
% LU.Y = column basis
% LU.IsLR = low-rank info

betatype = parameters.betatype;
blknum = parameters.numberOfBlocks;
variant = parameters.variant;

beta = zeros(blknum);
if betatype==1 || betatype==2
    % local threshold beta(i,j) = norm(A{i,j})
    for i=1:blknum
        for j=1:blknum
            beta(i,j) = mynorm(A{i,j});
        end
    end
elseif betatype==3  || betatype==4
    % global threshold beta(i,j) = norm(A)
    normA = mynorm(Afull);
    for i=1:blknum
        for j=1:blknum
            beta(i,j) = normA;
        end
    end
end

% initialize LU with A
LU = cell(blknum);
for i=1:blknum
    for j=1:blknum
        LU{i,j}.X = A{i,j};
        LU{i,j}.IsLR = 0;
        if strcmp(variant,'CUF')==1
            if i>j
                [LU{i,j}] = compress(LU{i,j}.X,beta(i,j),'L',parameters.lowRankThreshold);
            elseif i<j
                [LU{i,j}] = compress(LU{i,j}.X,beta(i,j),'U',parameters.lowRankThreshold);
            end
        end
    end
end

%%%%%% MAIN LOOP %%%%%%
waitBar = waitbar(0,'Running...');
for k=1:blknum
    waitbar(k/blknum,waitBar,['Low rank treshold 2 to ' num2str(log2(parameters.lowRankThreshold))]);
    %%%%%% COMPRESS step (UCF/CUF) %%%%%%
    if strcmp(variant,'UCF')==1
        for i=k+1:blknum
            [LU{i,k}] = compress(LU{i,k}.X,beta(i,k),'L',parameters.lowRankThreshold);
            [LU{k,i}] = compress(LU{k,i}.X,beta(k,i),'U',parameters.lowRankThreshold);
        end
    elseif strcmp(variant,'CUF')==1
        for i=k+1:blknum
            [LU{i,k}] = recompress(LU{i,k},beta(i,k),'L',parameters);
            [LU{k,i}] = recompress(LU{k,i},beta(k,i),'U',parameters);
        end
    end
    %%%%%% end COMPRESS step (UCF/CUF) %%%%%%

    %%%%%% FACTOR step %%%%%%
    LU{k,k}.X = mylu_nopiv(LU{k,k}.X);
    for i=k+1:blknum
        [LU{i,k}] = solve(LU{k,k}.X,LU{i,k},'L');
        [LU{k,i}] = solve(LU{k,k}.X,LU{k,i},'U');
    end
    %%%%%% end FACTOR step %%%%%%

    %%%%%% COMPRESS step (UFC/FUC) %%%%%%
    if strcmp(variant,'UFC')==1 || strcmp(variant,'FUC')==1
        if betatype==2 || betatype==4
            scalU = mynorm(triu(LU{k,k}.X));
            scalL = mynorm(tril(LU{k,k}.X,-1)+eye(length(LU{k,k}.X)));
        else
            scalU = 1; scalL = 1;
        end
    end
    if strcmp(variant,'UFC')==1
        for i=k+1:blknum
            [LU{i,k}] = compress(LU{i,k}.X,beta(i,k)/scalU,'L',parameters.lowRankThreshold);
            [LU{k,i}] = compress(LU{k,i}.X,beta(k,i)/scalL,'U',parameters.lowRankThreshold);
        end
    end
    %%%%%% end COMPRESS step (UFC/FUC) %%%%%%

    %%%%%% UPDATE step %%%%%%
    for i=k+1:blknum
        for j=k+1:blknum
            if i>=j
                [LU{i,j}] = update(LU{i,k},LU{k,j},LU{i,j},-1,'L',beta(i,j),parameters);
            else
                [LU{i,j}] = update(LU{i,k},LU{k,j},LU{i,j},-1,'U',beta(i,j),parameters);
            end
        end
    end
    %%%%%% end UPDATE step %%%%%%

    %%%%%% COMPRESS step (FUC) %%%%%%
    if strcmp(variant,'FUC')==1
        for i=k+1:blknum
            [LU{i,k}] = compress(LU{i,k}.X,beta(i,k)/scalU,'L',parameters.lowRankThreshold);
            [LU{k,i}] = compress(LU{k,i}.X,beta(k,i)/scalL,'U',parameters.lowRankThreshold);
        end
    end
    %%%%%% end COMPRESS step (FUC) %%%%%%
end
close(waitBar)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = solve(T,B,LorU)
% solves TZ = B for Z
Z = B;
L = tril(T,-1)+eye(length(T));
U = triu(T);
if B.IsLR && LorU=='L'
    Z.Y = B.Y/U;
    f = prod(size(U))*size(B.Y,1);
elseif LorU=='L'
    Z.X = B.X/U;
    f = prod(size(U))*size(B.X,1);
else
    Z.X = L\B.X;
    f = length(L)*(length(L)-1)*size(B.X,2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = update(A,B,C,sgn,LorU,betaM,parameters)
% computes Z = sgn*A*B + C
recomp = parameters.recomp;

% first, compute P=A*B
if A.IsLR && B.IsLR
    P.IsLR = 1;
    M = matmul(A.Y,B.X,parameters.precision);
    if recomp
        [M] = compress(M,betaM,LorU,parameters.lowRankThreshold,min(size(M))-1);
    else
        Mtemp.X = M;
        Mtemp.IsLR = 0;
        M = Mtemp;
    end
    if M.IsLR
        P.X = matmul(A.X,M.X,parameters.precision); P.Y = matmul(M.Y,B.Y,parameters.precision);
    else
        if size(M,1) < size(M,2)
            P.X = A.X; P.Y = matmul(M.X,B.Y,parameters.precision);
        else
            P.X = matmul(A.X,M.X,parameters.precision); P.Y = B.Y;
        end
    end
elseif A.IsLR
    P.IsLR = 1;
    P.X = A.X; P.Y = matmul(A.Y,B.X,parameters.precision);
elseif B.IsLR
    P.IsLR = 1;
    P.X = matmul(A.X,B.X,parameters.precision); P.Y = B.Y;
else
    P.IsLR = 0;
    P.X = matmul(A.X,B.X,parameters.precision);
end

% then, compute Z=sgn*P+C
if C.IsLR && P.IsLR
    Z.X = [C.X,sgn*P.X]; Z.Y = [C.Y;P.Y];
    Z.IsLR = 1;
elseif C.IsLR
    Z.X = matmul(C.X,C.Y,parameters.precision)+sgn*P.X;
    Z.IsLR = 0;
elseif P.IsLR
    Z.X = C.X+sgn*matmul(P.X,P.Y,parameters.precision);
    Z.IsLR = 0;
else
    Z.X = C.X+sgn*P.X;
    Z.IsLR = 0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = recompress(A,beta,LorU,parameters)
if ~A.IsLR || size(A.X,2)==0
    B = A;
    return
end
B.IsLR = 1;
maxrk = size(A.X,1)*size(A.Y,2)/(size(A.X,1)+size(A.Y,2));
if LorU=='U'
    C = compress(A.X,beta,'U',parameters.lowRankThreshold,min(size(A.X)));
    B.X = C.X; B.Y = matmul(C.Y,A.Y,parameters.precision);
    %    C = compress(A.X,beta,'U',min(size(A.X)));
    %    D = matmul(C.Y,A.Y);
    %    E = compress(D,beta,'L',min(size(D)));
    %    B.X = matmul(C.X,E.X); B.Y = E.Y;
else
    C = compress(A.Y,beta,'L',parameters.lowRankThreshold,min(size(A.Y)));
    B.X = matmul(A.X,C.X,parameters.precision); B.Y = C.Y;
    %    C = compress(A.Y,beta,'L',min(size(A.Y)));
    %    D = matmul(A.X,C.X);
    %    E = compress(D,beta,'L',min(size(D)));
    %    B.X = E.X; B.Y = matmul(E.Y,C.Y);
end
if size(B.X,2) > maxrk
    B.IsLR = 0;
    B.X = matmul(B.X,B.Y);
    B.Y = [];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = matmul(A,B,uprec)
if uprec=='h'
    C = matmul_half(A,B);
else
    C = A*B;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = matmul_half(A,B)
C = zeros(size(A,1),size(B,2));
for k=1:size(A,2)
    C = chop(C + chop(A(:,k)*B(k,:), 11), 11);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = compress(A,beta,LorU,lrth,kmax)
% compress based on truncated SVD
if nargin<5
    s=size(A);
    kmax=prod(s)/sum(s);
end
[U,Smat,V] = svd(full(A));
Svec = diag(Smat);
if isa(U,'double')
    U = sparse(U);
    Smat = sparse(Smat);
    V = sparse(V);
end

rk = 1;
n = length(Svec);
while rk<min(n,kmax+1) && sqrt(sum(Svec(rk+1:end).^2))>lrth*beta
    rk = rk+1;
end
if (rk <= kmax)
    B.IsLR = 1;
    if LorU=='L'
        B.X = U(:,1:rk);
        B.Y = Smat(1:rk,1:rk)*V(:,1:rk)';
    else
        B.X = U(:,1:rk)*Smat(1:rk,1:rk);
        B.Y = V(:,1:rk)';
    end
else
    B.IsLR = 0;
    B.X = A;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = mylu_nopiv(A)
% lu with no pivoting
n = size(A, 1);
if isa(A,'double')
    L = speye(n);
else
    L = eye(n);
end
for k = 1 : n
    L(k + 1 : n, k) = A(k + 1 : n, k) / A(k, k);
    for l = k + 1 : n
        A(l, :) = A(l, :) - L(l, k) * A(k, :);
    end
end
U = A;
B = tril(L,-1) + U;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function my_fprintf(fout,args)
% printf to either file or terminal
if fout == 0
    fprintf(args{1:end});
else
    fprintf(fout,args{1:end});
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = bwderr(Afull,S,blknum)
xtrue = ones(size(Afull,1),1);
rhs = Afull*xtrue;

y = fwd_subst(S,blknum,rhs);
x = bwd_subst(S,blknum,y);
normA = mynorm(Afull);
normx = mynorm(xtrue);
err = mynorm(Afull*x-rhs)/(normA*normx);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fwd_subst(S,blknum,rhs)
y=rhs;
jend = 0;
for j=1:blknum
    jbeg = jend+1;
    jend = jbeg + size(S{j,j}.X,1) -1;
    Ljj = tril(S{j,j}.X,-1)+speye(size(S{j,j}.X));
    y(jbeg:jend,:) = Ljj\y(jbeg:jend,:);
    iend = jend;
    for i=j+1:blknum
        ibeg = iend+1;
        iend = ibeg + size(S{i,j}.X,1) -1;
        if S{i,j}.IsLR
            y(ibeg:iend,:) = y(ibeg:iend,:) - S{i,j}.X*(S{i,j}.Y*y(jbeg:jend,:));
        else
            y(ibeg:iend,:) = y(ibeg:iend,:) - S{i,j}.X*y(jbeg:jend,:);
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = bwd_subst(S,blknum,y)
x=y;
ibeg = size(y,1)+1;
for i=blknum:-1:1
    iend = ibeg-1;
    ibeg = iend - size(S{i,i}.X,1) +1;
    jend = iend;
    for j=i+1:blknum
        jbeg = jend+1;
        if S{i,j}.IsLR
            jend = jbeg + size(S{i,j}.Y,2) -1;
            x(ibeg:iend,:) = x(ibeg:iend,:) - S{i,j}.X*(S{i,j}.Y*x(jbeg:jend,:));
        else
            jend = jbeg + size(S{i,j}.X,2) -1;
            x(ibeg:iend,:) = x(ibeg:iend,:) - S{i,j}.X*x(jbeg:jend,:);
        end
    end
    Uii = triu(S{i,i}.X);
    x(ibeg:iend,:) = Uii\x(ibeg:iend,:);
end
end



