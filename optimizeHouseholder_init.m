function theVs = optimizeHouseholder_init(Q, hh)
[n, n] = size(Q);

theVs = [];
minvs2 = [];
[Vf, Df] = eig(Q);
df_diag = diag(Df);

visited = [];
for index = 1:n
    if (ismember(index, visited))
        continue;
    end
    val = df_diag(index);
    indices = find( abs(df_diag-val) < 10e-7 );
    
    if (length(indices) > 1)
        W = Vf(:, indices);
        [W, R] = qr(W);
        W = W(:, 1:length(indices));
        Vf(:, indices) = W;
    end
    visited = union(visited, indices);
end
Df = Vf'*Q*Vf;

Qnew = Q;
while (size(theVs,2) < hh)
    [min_v, min_i] = min(diag(real(Df)));
    if (min_v >0)
        break;
    end

    if (abs(Df(min_i, min_i) - real(Df(min_i, min_i))) < 10e-7)
        theVs = [theVs Vf(:, min_i)];
        minvs2 = [minvs2 real(Df(min_i, min_i))];
        Df(min_i, min_i) = 2;
        
    else
        v1 = Vf(:, min_i);
        v2 = Vf(:, min_i+1);
        u1 = 1/sqrt(2)*(v1+v2);
        U1 = eye(n) - 2*(u1*u1');
        Qnew = Qnew*U1;
        
        theVs = [theVs u1];
        minvs2 = [minvs2 2*real(Df(min_i, min_i))];
        
        lmin = Df(min_i, min_i)+Df(min_i+1, min_i+1);
        
        g = -sqrt(lmin/2+1)/2;
        f = -sqrt(-lmin/2+1)/2;
        
        u2 = (g+1i*f)*v1 + (g-1i*f)*v2;
        U2 = eye(n) - 2*(u2*u2');
        Qnew = Qnew*U2;
        
        theVs = [theVs u2];
        minvs2 = [minvs2 -2];
        
        Df(min_i, min_i) = 2;
        Df(min_i+1, min_i+1) = 2;
    end
end

T = eye(n);
for ii = 1:size(theVs,2)
%     T = (eye(n) - 2*theVs(:, ii)*theVs(:, ii)')*T;
    T = apply_reflector_left(theVs(:, ii), T);
end

theVs = theVs(:, 1:hh);
