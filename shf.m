function [reflectors, s, d, val, U] = shf(S, h, changeSpectrum, changeD, s, reflectors)
%% The Symmetric Householder Factorization algorithm for symmetric matrices, Algorithm 1 in the paper

[n, ~] = size(S);

D = eye(n); d = diag(D);

K = 150;
val = inf;
total = norm(S,'fro')^2;
for i = 1:K
    Ak = D*S*D;
    for j = 1:h
%         Ak = (eye(n) - 2*reflectors(:,j)*reflectors(:,j)')*Ak*(eye(n) - 2*reflectors(:,j)*reflectors(:,j)');
        Ak = apply_reflector_left(reflectors(:,j), Ak);
        Ak = apply_reflector_right(reflectors(:,j), Ak);
    end
    Bk = diag(s);
    
    % update reflectors
    U = eye(n);
    for k = h:-1:1
%         Ak = (eye(n) - 2*reflectors(:,k)*reflectors(:,k)')*Ak*(eye(n) - 2*reflectors(:,k)*reflectors(:,k)');
        Ak = apply_reflector_left(reflectors(:,k), Ak);
        Ak = apply_reflector_right(reflectors(:,k), Ak);
        
        reflectors(:,k) = get_uk(Ak, Bk, reflectors(:,k));
        
%         U = (eye(n) - 2*reflectors(:,k)*reflectors(:,k)')*U;
        U = apply_reflector_left(reflectors(:,k), U);
        
%         Bk = (eye(n) - 2*reflectors(:,k)*reflectors(:,k)')*Bk*(eye(n) - 2*reflectors(:,k)*reflectors(:,k)');
        Bk = apply_reflector_left(reflectors(:,k), Bk);
        Bk = apply_reflector_right(reflectors(:,k), Bk);
    end
    
    % update s
    if (changeSpectrum)
        s = diag(U'*S*U);
        Bk = U*diag(s)*U';
    end
    
    % update D
    if (i > 1 && changeD)
        for k = 1:n
            stilde = S(:,k);
            stilde(k) = [];

            btilde = Bk(:,k);
            btilde(k) = [];

            if (norm(stilde - btilde) < norm(stilde + btilde))
                D(k,k) = -1;
            else
                D(k,k) = 1;
            end
        end
        d = diag(D);
    end
    old_val = val;
    val = (norm(S-D*Bk*D,'fro')^2)/total;
    
    if (abs(val-old_val) < 10e-8)
        break;
    end
end

