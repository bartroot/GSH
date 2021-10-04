function [B] = Europe_centered(A)
%
% centers Europe in the matrix 

B = zeros(size(A));
B(:,1:end/2) = A(:,end/2+1:end);
B(:,end/2+1:end) = A(:,1:end/2);