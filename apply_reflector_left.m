function y = apply_reflector_left(reflector, X)
% y = (eye(n) - 2*reflector*reflector')*X;
y = reflector'*X; y = reflector*y;
y = X - 2*y;
