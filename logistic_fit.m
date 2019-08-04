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

function [P, E] = logistic_fit(F, t, k, LB, UB, Dfit, Nfit)
% Fit every row in the matrix to a logistic function
% 
% F = an expression matrix [genes]x[time]]
% t = time points (in min)
% k = optimization repeats
% LB,UB = parameter bounds
%
% Output:
%  P = <b,h0,h1,t1> (parameters)
%  E = Error

if (nargin < 3)
    k = 100;
end
if (nargin < 5)
    LB = [0 0 0 0];
    UB = [10 inf inf max(t)];
end
if (nargin < 6)
    Dfit = 1e-5;
end
if (nargin < 7)
    Nfit = 1e3;
end

OPTS = optimset('TolX',Dfit, 'TolFun',Dfit, ...
    'MaxFunEvals',Nfit, 'MaxIter', Nfit, ...
    'Jacobian','on', 'Display','off');

n = size(F,1);
P = zeros(n, 4);
E = zeros(n, 1);
for i = 1:n
    Pi = zeros(1, 4);
    Ei = inf;
    P0 = logistic_rnd_init(F(i,:), t, k, LB, UB);
    for j = 1:k
        [p, e] = lsqnonlin(@errf, P0(j,:), LB, UB, OPTS, t, F(i,:));
        if (e < Ei)
            Ei = e;
            Pi = p;
        end
    end
    P(i,:) = Pi;
    E(i,:) = Ei;
end


function [e,d] = errf(P, t, F)
    
f = logistic_eval(P, t);
e = (f - F);
if (nargout > 1)
    d = logistic_df(P, t);
    d = d';
end
function d = logistic_df(P, x, s, p)
% Logistic function gradient with respect to all paramenters,
%  with adjusted size
%
% P = parameter vector: <b,h0,h1,t1>
% x = time points
% s = number of rows in final gradient vector (d)
% p = start position of the logistic gradient within the vector
%
% d = <d/db, d/dh0, d/dh1, d/dt1>
%
% Logistic function:
%  f(x) = h0 + (h1 - h0)/(1 + exp( b*(x - t1)));

if (nargin < 3)
    s = 4;
end
if (nargin < 4)
    p = 1;
end


% estimate gradient
gr = logistic_grad(P, x);

% convert vector to the correct size
[n,m] = size(gr);
d = zeros(s,m);
d(p:(p+n-1),:) = gr;


function d = logistic_grad(P, x)

% parameters
b = P(:,1);
h0 = P(:,2);
h1 = P(:,3);
t1 = P(:,4);
s = sigmoid(b, t1, x);

% derivative
%   S = 1/(1+E)  ==>  E/(1+E)^2 = E*(S^2) = (1-S)*S
Db = (1 - s) .* s .* (h1 - h0) .* (t1 - x);
Dh0 = 1 - s;
Dh1 = s;
Dt1 = (1 - s) .* s .* (h1 - h0) .* b;

% final
d = [Db; Dh0; Dh1; Dt1];

function p0 = logistic_rnd_init(F, t, k, LB, UB)
% init random parameters 
% 4 parameters: <b,h0,h1,t1>
% produce <k> inits per line in F

if (nargin < 3)
    k = 1;
end
if (nargin < 5)
    LB = (-1)*inf(1,4);
    UB = inf(1,4);
end

mnt = max(min(t), LB(4));
mxt = min(max(t), UB(4));

n = size(F,1);
p0 = zeros(n*k,4);
for i=1:n
    % t1
    t1 = mnt(1) + (mxt(1) - mnt(1))*rand(k,1);

    % h0,h1
    mnf = min(F(i,:));
    mxf = max(F(i,:));
    r = mxf - mnf;
    mnf = max([LB(2:3); repmat(mnf - 0.5*r,1,2)]);
    mxf = min([UB(2:3); repmat(mxf + 0.5*r,1,2)]);
    h0 = mnf(1) + (mxf(1) - mnf(1))*rand(k,1);
    h1 = mnf(2) + (mxf(2) - mnf(2))*rand(k,1);
    
    % b
    slp = abs((F(i,2:end) - F(i,1:end-1))./(t(2:end) - t(1:end-1))); % slope
    abs_b = 4*max(slp)/r;
    if (r == 0)
        abs_b = 0;
    end
    mnb = max(-2*abs_b, LB(1));
    mxb = min(2*abs_b, UB(1));
    b = mnb + (mxb - mnb)*rand(k,1);
    
    p0((i-1)*k+1:i*k,:) = [b,h0,h1,t1];
end

p0 = unique(p0,'rows');
