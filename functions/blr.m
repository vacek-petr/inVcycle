function [approx,info] = blr(A,b,BLRLU)
% BLR direct solver based on BLR LU decompostion (see e.g. Higham, Mary 2021)
% parts of code taken from
% https://gitlab.com/theo.andreas.mary/BLRstability
%
% Inputs:
% A - matrix
% b - right-hand side vector
% BLRLU - precomputed BLR LU decomposition of A by the script
% computeBLRLU.m
%
% Outputs:
% approx - computed approximation of A^(-1)*b
% info - structure - contains
%   res2norm - Euclidian norm of residual b-A*approx
%   relres2norm - Euclidian norm of residual b-A*approx scaled by Euclidian
%   norm of b
%   numberOfIterations = 1


numberOfBlocks = size(BLRLU,1);

y = fwd_subst(BLRLU,numberOfBlocks,b);
approx = bwd_subst(BLRLU,numberOfBlocks,y);

info.res2norm = norm(b-A*approx);
info.relres2norm = info.res2norm/norm(b);
info.approx2norm = norm(approx);
info.rhs2norm = norm(b);
info.numberOfIterations = 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
