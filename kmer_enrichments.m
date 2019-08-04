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

function P = kmer_enrichments(I, S, tail, test, minI)
% Input:
%  I = binary group assignments (kmers x UTRs)
%  S = scores to test, NaN values are ignored (UTRs x 1)
%  tail = 0: both, 1: right (mean > BG), -1: left (mean < BG)
%  test = 0: t-test
%         1: ks-test
%         2: permutation test (SLOW)
%         3: hypergeometric pvalue (discrete S)
%        -1: mutual information
%
% Output:
%  P = pvalue

if (nargin < 3)
    tail = 0;
end
if (nargin < 4)
    test = 0;
end
if (nargin < 5)
    minI = 10;
end

% discrete S
if (test == 3)
    u = unique(S(~isnan(S)));
    n = max(size(u));
    J = zeros(size(S,1),n);
    for k = 1:n
        J(:,k) = (S==u(k));
    end
    J(isnan(S),:) = NaN;
    S = J;
end
n = size(S,2);
m = size(I,1);

% pvalues
P = ones(m,n);
for k = 1:n
    j = (~isnan(S(:,k)));
    i = sum(I(:,j),2)>1;
    if (test >= 0)
        P(:,k) = ones(m,1);
        P(i,k) = pvdist(I(i,j), [], S(j,k), tail, test, minI);
    else
        P(:,k) = zeros(m,1);
        P(i,k) = midist(I(i,j), [], S(j,k), minI);
    end
end


function d = pvdist(x, y, S, tail, test, minI)
% p-values for comparing the following hypothesis:
%   H0 = all data in S comes from the same distribution
%   H1 = data in S marked by "1" (in either x or y) is different
%        distribution than the rest of S
%
%  tail = 0: both, 1: right (mean > BG), -1: left (mean < BG)
%  test = 0: t-test, 1: permutation test (SLOW), 2: ks-test, 3: hypergeometric

if (isempty(y))
    y = zeros(1,size(x,2));
end

if ((size(x,2) ~= size(S,1)) || (size(y,2) ~= size(S,1)))
    d = NaN;
    return;
end

M = size(S,1);
N = sum(S>0);

nx = size(x,1);
ny = size(y,1);
d = inf(nx,ny);
for ix = 1:nx
    for iy = 1:ny
        J = x(ix,:)' + y(iy,:)';
        bg = S(J>=0);
        dt = S(J> 0);
        
        if (min(size(bg,1),size(dt,1)) < minI)
            d(ix,iy) = 1;
            continue;
        end
        
        if (test == 3) % hypergeometric pvalue
            X = sum(dt>0);
            K = sum(J);
            d(ix,iy) = hypergeometric_pvalue(X,M,K,N);
        
        elseif (test == 2) % permutation test
            [~,d(ix,iy)] = permtest(dt, bg, 1000, tail);
            
        elseif (test == 1) % ks-test
            if (tail > 0)
                [~,d(ix,iy)] = kstest2(dt, bg, 'Tail','smaller'); % dt > BG
            elseif (tail < 0)
                [~,d(ix,iy)] = kstest2(dt, bg, 'Tail','larger'); % dt < BG
            else
                [~,d(ix,iy)] = kstest2(dt, bg);
            end
            
        else % t-test
            z = (dt-mean(bg))./std(bg);
            if (tail > 0)
                [~,d(ix,iy)] = ttest(z,[],'tail','right'); % mean > BG
            elseif (tail < 0)
                [~,d(ix,iy)] = ttest(z,[],'tail','left'); % mean < BG
            else
                [~,d(ix,iy)] = ttest(z);
            end
        end
    end
end


function d = midist(x, y, S, minI)
% mutual information between rows of x to columns of S

if (isempty(y))
    y = zeros(1,size(x,2));
end

if ((size(x,2) ~= size(S,1)) || (size(y,2) ~= size(S,1)))
    d = NaN;
    return;
end

mxS = max(S(:));
mnS = min(S(:));
nBin = 100;
discreteS = mxS * round(nBin*(S - mnS)./mxS)/nBin + mnS;

nx = size(x,1);
ny = size(y,1);
d = inf(nx,ny);
for ix = 1:nx
    for iy = 1:ny
        J = x(ix,:)' + y(iy,:)';
        bg = S(J>=0);
        dt = S(J> 0);
        if (min(size(bg,1),size(dt,1)) < minI)
            d(ix,iy) = inf;
            continue;
        end
        
        d(ix,iy) = mutual_information(J,discreteS);
    end
end

function [h,p] = permtest(X,Y,N,tail)
% test the null hypothesis that X was sampled from a distribution
% with a particular mean
% more general than parametric t-tests in that it does not assume that
% the data were sampled from a normal distribution
%
% X = test set
% Y = data sampled from the bacground distribution
% N = number of permutations

if (nargin < 3)
    N = 1000;
end
if (nargin < 4)
    tail = 0;
end

% test statistics
D = [X;Y];
h = mean(X) - mean(D);

% use random permutations to calculate an empirical pvalue
k = max(size(X));
M = mean_diff(D,N,k);
if (tail > 0)
    p = (1+10*sum(M>h))./(1+10*N);
elseif (tail < 0) 
    p = (1+10*sum(M<h))./(1+10*N);
else
    p = (1+10*sum(abs(M)>abs(h)))./(1+10*N);
end


function M = mean_diff(D,N,k)
% difference of means of N random samples of size K

rng('shuffle');
n = max(size(D));

M = zeros(N,1);
m = mean(D);
for i = 1:N
    p = randperm(n,k);
    M(i) = mean(D(p)) - m;
end

function M = mutual_information(A,B)
% discrete mutual information (log2 base)

% joint probability
[uAB,~,t] = unique([A B],'rows');
nAB = accumarray(t,1);
pAB = nAB./sum(nAB);

% marginalize A's probability
[uA,~,i] = unique(uAB(:,1));
pA = accumarray(i,pAB);
pA = pA(i);

% marginalize B's probability
[uB,~,j] = unique(uAB(:,2));
pB = accumarray(j,pAB);
pB = pB(j);

% mutual information = sum_x sum_y p(x,y)*log(p(x,y)/p(x)p(y))
logP = log2(pAB) - log2(pA) - log2(pB);
M = sum(pAB.*logP);
if (M<0)
    M = 0;
end

function p = hypergeometric_pvalue(X,M,K,N,lower)
% p = hypergeometric_pvalue(X,M,K,N,lower)
%
% X = number of items with the desired characteristic in sample
% M = size of the population
% K = number of items with the desired characteristic in population
% N = size of sample drawn
% lower = calculate the lower tail (by default calculates the upper tail)
%
% Upper tail = p is the probability of drawing X or MORE items with the 
%              desired characteristic when drawing N out of M items in 
%              which K items has the desired characteristic
% Lower tail = p is the probability of drawing X or LESS items with the 
%              desired characteristic when drawing N out of M items in 
%              which K items has the desired characteristic

if (nargin < 5)
    lower = 0;
end

if (lower) 
    p = hygecdf(X,M,K,N);
else
    p = hygecdf(X,M,K,N,'upper') + hygepdf(X,M,K,N);
end
