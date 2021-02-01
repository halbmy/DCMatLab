function [U,sm,X,V] = gensvd(A,C) 

% GENSVD - Generalized SVD of matrix pair A,C
% calls routine gsvd (matlab built-in)
% works for both size(A,1) greater/less than size(A,2)
% [U,sm,X,V] = gensvd(A,C)

[d,m]=size(A);
if d>=m, % regular procedure
    [U,V,W,s,c] = gsvd(A,full(C),0); 
else % switching A and C for svd
    [V,U,W,c,s] = gsvd(full(C),A,0); 
end
p=min(size(A));
% singular values
sm = [diag(s(1:p,1:p)),diag(c(1:p,1:p))]; 

if (nargout < 2), % only generalized singular values
    U = sm(:,1)./sm(:,2); 
else % complete matrices
    X = inv(W'); 
end