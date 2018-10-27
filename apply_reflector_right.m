function y = apply_reflector_right(reflector, X)
% y = X*(eye(n) - 2*reflector*reflector');
y = X*reflector; y = y*reflector';
y = X - 2*y;
