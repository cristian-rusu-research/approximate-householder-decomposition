function [u, branch1, branch2, iters, vals] = get_uk(A, B, u)
[n, ~] = size(A);

W = A*B+B*A;
W = (W+W')/2;
A = (A+A')/2;
B = (B+B')/2;
if (nargin == 2 || norm(u) <= 10e-10)
    [u1, ~, ~] = eigs(W, 1, 'SA');

    [vb_plus, db_plus, ~] = eigs(B, 1, 'LA'); [vb_minus, db_minus, ~] = eigs(B, 1, 'SA');
    [va_plus, da_plus, ~] = eigs(A, 1, 'LA'); [va_minus, da_minus, ~] = eigs(A, 1, 'SA');

    vals = [da_plus*db_plus da_minus*db_plus da_plus*db_minus da_minus*db_minus];
    [maxv, maxi] = max(vals);
    if (maxi == 1)
        vv = kron(vb_plus, va_plus);
    end
    if (maxi == 2)
        vv = kron(vb_plus, va_minus);
    end
    if (maxi == 3)
        vv = kron(vb_minus, va_plus);
    end
    if (maxi == 4)
        vv = kron(vb_minus, va_minus);
    end

    V = reshape(vv', [n n]);
    V = V + V';
    [U1, ~, ~] = svds(V, 1);
    u2 = U1(:,1);
    
    u2 = u2 - u1'*u2*u1; u2 = u2/norm(u2);

    Au1 = A*u1; Au2 = A*u2;
    Bu1 = B*u1; Bu2 = B*u2;
    a = u1'*Au1; b = u2'*Au2; c = u1'*Au2;
    x = u1'*Bu1; y = u2'*Bu2; z = u1'*Bu2;
    Wu1 = W*u1; Wu2 = W*u2;
    h = u1'*Wu1; j = u2'*Wu2; k = u1'*Wu2;

    fun1 = @(e)(((1-e^2/2)^2*h + ((e^2-e^4/4)^(1/2))^2*j + 2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*k) -2*((1-e^2/2)^2*a+((e^2-e^4/4)^(1/2))^2*b+2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*c)*((1-e^2/2)^2*x+((e^2-e^4/4)^(1/2))^2*y+2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*z));
    fun2 = @(e)(((1-e^2/2)^2*h + ((e^2-e^4/4)^(1/2))^2*j - 2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*k) -2*((1-e^2/2)^2*a+((e^2-e^4/4)^(1/2))^2*b-2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*c)*((1-e^2/2)^2*x+((e^2-e^4/4)^(1/2))^2*y-2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*z));
    [X1, FVAL1, EXITFLAG1] = fminbnd(fun1, 0, sqrt(2));
    [X2, FVAL2, EXITFLAG2] = fminbnd(fun2, 0, sqrt(2));
    if (FVAL1 < FVAL2)
        epsilon = X1;
        u = (1-epsilon^2/2)*u1 + (epsilon^2-epsilon^4/4)^(1/2)*u2;
    else
        epsilon = X2;
        u = (1-epsilon^2/2)*u1 - (epsilon^2-epsilon^4/4)^(1/2)*u2;
    end
end

vals = u'*W*u - 2*(u'*A*u)*(u'*B*u);
iters = 0;
branch1 = 0; branch2 = 0;
while(1)
    Au = A*u; Bu = B*u;
    g = 2*W*u - 4*u'*Au*Bu - 4*u'*Bu*Au;
    g = g - (u'*g)*u;
    g = g/norm(g);

    Ag = A*g;
    Bg = B*g;
    a = u'*Au; b = g'*Ag; c = u'*Ag;
    x = u'*Bu; y = g'*Bg; z = u'*Bg;
    Wu = W*u; Wg = W*g;
    h = u'*Wu; j = g'*Wg; k = u'*Wg;

    fun1 = @(e)(((1-e^2/2)^2*h + (e^2-e^4/4)*j + 2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*k) -2*((1-e^2/2)^2*a+(e^2-e^4/4)*b+2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*c)*((1-e^2/2)^2*x+(e^2-e^4/4)*y+2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*z));
    fun2 = @(e)(((1-e^2/2)^2*h + (e^2-e^4/4)*j - 2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*k) -2*((1-e^2/2)^2*a+(e^2-e^4/4)*b-2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*c)*((1-e^2/2)^2*x+(e^2-e^4/4)*y-2*(1-e^2/2)*((e^2-e^4/4)^(1/2))*z));
    [X1, FVAL1, EXITFLAG1] = fminbnd(fun1, 0, sqrt(2));
    [X2, FVAL2, EXITFLAG2] = fminbnd(fun2, 0, sqrt(2));
    old_val = u;
    
    if (vals(end)<min(FVAL1, FVAL2))
        break;
    end
    
    if (FVAL1 < FVAL2)
        epsilon = X1;
        u = (1-(epsilon^2)/2)*u + (epsilon^2-(epsilon^4)/4)^(1/2)*g;
        branch1 = branch1 + 1;
    else
        epsilon = X2;
        u = (1-(epsilon^2)/2)*u - (epsilon^2-(epsilon^4)/4)^(1/2)*g;
        branch2 = branch2 + 1;
    end

    val = u'*W*u - 2*(u'*A*u)*(u'*B*u);
    if (vals(end) < val)
        u = old_val;
        break;
    end
    if ( abs(norm(u)-1) > 10e-7)
        u = old_val;
        break;
    end
    vals = [vals val];
    
    iters = iters + 1;
    
    if (norm(old_val-u)^2 <= 10e-5)
        break;
    end
end
