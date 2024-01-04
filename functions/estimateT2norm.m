function T2norm = estimateT2norm(mh)
% T2NORM function for estimating the Euclidean norm of the matrix T

A=mh.A;  P=mh.P; J = mh.numberOfLevels;


T2norm = sqrt(eigs(@(x) T_transpose_apply(T_apply(x,mh),mh), size(mh.A{mh.numberOfLevels},2),1,'largestabs',Tolerance=1e-6));

function x = T_apply(x,mh)
    A=mh.A; P=mh.P; J=mh.numberOfLevels;

    % finest level smoothing
    U = triu(A{J});
    L = tril(A{J});
    D = diag(diag(A{J}));
    x = x - A{J}*(U\(D*(L\x)));

    for j=J-1:-1:2
        x = P{j+1}'*x;
        U = triu(A{j});
        L = tril(A{j});
        D = diag(diag(A{j}));
        x = x - A{j}*(U\(D*(L\x)));
    end
    x = P{2}'*x;
end

function x = T_transpose_apply(x,mh)
    A=mh.A; P=mh.P; J=mh.numberOfLevels;

    for j=2:J-1
        x = P{j}*x;
        U = triu(A{j});
        L = tril(A{j});
        D = diag(diag(A{j}));

        x = x - (U\(D*(L\(A{j}*x))));
    end
    x = P{J}*x;
    % finest level smoothing
    U = triu(A{J});
    L = tril(A{J});
    D = diag(diag(A{J}));
    x = x - U\(D*(L\(A{J}*x)));

end



end
