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

function [W,S,b0] = kmer_positional_weights(ids,seq,model_file)
% calculate positional weights by kmer linear regression model
%
% Input:
%  ids = sequence ids
%  model_file = kmer linear regression model
%
% Output:
%  W[i] = sum of weights of all kmers that overlap the i-th position in sequence
%        (rows = ids, columns = position on sequence)
%     S = overall sum of all kmer weights


% model data
load(model_file,'Krows','B0','b0');
Brows = Krows;
Blen = cellfun(@length,Brows);

% weights
m = size(ids,1);
n = size(seq,2);

W = zeros(m,n);
S = zeros(m,1);
u = unique(Blen);
u = u(u<=n);
for i = 1:m
    S(i) = b0;
    for k = 1:max(size(u))
        [wk,sk] = calculate_weights(seq(i,:),Brows(Blen==u(k)),B0(Blen==u(k)),u(k));
        W(i,:) = W(i,:) + wk;
        S(i) = S(i) + sk;
    end
end



function [W,S] = calculate_weights(seq,Brows,Bweights,k)
% W[i] = sum of all kmers that overlap the i-th position in sequence
% S = sum of all kmer weights

n = length(seq);
W = zeros(1,n);

% split sequence into kmers of size k
K = strings(n,1);
T = repmat(char(0),k,n-k+1);
for i = 1:k
    T(i,:) = seq(1+(i-1):n-(k-i));
end
K(1:n-k+1) = cellstr(T');

% kmer weights
KW = zeros(1,n);
[i,j] = ismember(K,Brows);
KW(i) = Bweights(j(j>0));
S = sum(KW);

% sum kmers per position
for i = 1:k
    W = W + [zeros(1,i-1) KW(1:n-i+1)];
end
