% ------------------------------------------------------------------------
% Copyright 2019 Michal Rabani
% 
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, 
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice, 
%    this list of conditions and the following disclaimer in the documentation 
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ------------------------------------------------------------------------

function [A,B,E] = linear_fit(X, Y, solvertype)
% solve a matrix linear regression problem
% Y = X*A + B

if (nargin < 3)
    solvertype = 0;
end

% add an intercept
if (nargout > 1)
    n = size(X,1);
    Z = [ones(n,1) X];
else
    Z = X;
end

if (solvertype == 1)
    % least squares optimal solution
    A = (Z'*Z)\(Z'*Y);
elseif (solvertype == 2)
    % can handle rank deficient matrices
    [C,R] = qr(sparse(Z),Y,'vector');
    A = R\C;
elseif (solvertype == 3)
    A = regress(Y,Z);
elseif (solvertype == 4)
    % (specific case of mldivide)
    A = linsolve(Z,Y);
elseif (solvertype == 5)
    % lasso regularization
    [a,fitinfo] = lasso(Z,Y,'CV',10);
    i = fitinfo.Index1SE;
    A = a(:,i);
else
    A = mldivide(Z,Y);
end

% intercept
if (nargout > 1)
    B = A(1);
    A = A(2:end);
else
    B = 0;
end

% least squares error
y = X*A + B;
E = sum((Y-y).^2);
