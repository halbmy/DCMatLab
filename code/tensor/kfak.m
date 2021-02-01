function k=kfak(A,B,M,N)

% k = kfak(A,B,M,N)

k=2*pi./(1/norm(A-M)-1/norm(A-N)-1/norm(B-M)+1/norm(B-N));