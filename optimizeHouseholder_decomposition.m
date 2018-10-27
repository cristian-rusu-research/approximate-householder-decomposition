function [U, X1, X2, theVs, tus, err, theVsoriginal] = optimizeHouseholder_decomposition(Q, h)
%% Householder decomposition of an orthonormal matrix, Section II of the paper

tic;
[n, ~] = size(Q);

% an initialization of Householder reflectors by QR decomposition
% [AA, theUs] = QR_tria_hh(Data, h);
% D = eye(n);
% for ii = 1:n
%     if (AA(ii,ii) < 0)
%         D(ii,ii) = -1;
%     end
% end
% theUs = fliplr(theUs);

theVs = optimizeHouseholder_init(Q, h);

theVsoriginal = theVs;
U = eye(n);
for j = 1:h
%     U = (eye(n) - 2*(theUs(:,j)*theUs(:,j)'))*U;
    U = apply_reflector_left(theVs(:,j), U);
end

X1 = eye(n);
X2 = eye(n);
err = norm(Q-X1*U*X2, 'fro')^2/norm(Q, 'fro')^2*100;

UU = X1*U*X2;
for iii = 1:n
    if norm(Q(:,iii)-UU(:,iii)) > norm(Q(:,iii)+UU(:,iii))
        X2(iii,iii) = -X2(iii,iii);
    end
end
err = [err norm(Q-X1*U*X2, 'fro')^2/norm(Q, 'fro')^2*100];

tus = toc;
