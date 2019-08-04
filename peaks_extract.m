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

function [X,Xu,Xw,Y,Yu,Yw,Z,Zu,Zw] = peaks_extract(id,seq,w,min_peak_weight,min_peak_width,peak_extend)
% find peaks in input weights
%
% Output:
%  X/Y = all peaks: <peak weight> <reporter id> <start> <end> <sequence>
%        peak weight = sum all weights of peak positions
%  Xu/Yu = unique peaks: <peak weight> <sequence>
%  Xw/Yw = unique peak sequences: <peak weight> <sequence>
%
% X,Xu,Xw = positive peaks
% Y,Yu,Yw = negative peaks

if (nargin < 5)
    min_peak_width = 1;
end
if (nargin < 6)
    peak_extend = 2;
end

% positive peaks
X = select_peaks(id,w,seq,1,min_peak_weight,min_peak_width,peak_extend);
np = size(X,1);
fprintf('  positive peaks: %d peaks identified \n', np);
if (np > 0)
    Xu = unique_peaks(X(:,[1 5]));
    fprintf('  positive peaks: %d unique peaks identified (%.1f%%)\n',size(Xu,1),100*size(Xu,1)./np);
    plen = cellfun(@length,Xu(:,2));
    pscore = abs(cell2mat(Xu(:,1))./plen);
    fprintf('    mean length = %.1f (min=%d, max=%d)\n', mean(plen),min(plen),max(plen));
    fprintf('    mean absolute height = %.2f (min=%.2f, max=%.2f)\n', mean(pscore),min(pscore),max(pscore));
    Xw = unique_peak_sequences(Xu);
    fprintf('  positive peaks: %d unique peak sequences (%.1f%%)\n\n', size(Xw,1),100*size(Xw,1)./np);
else
    Xu = [];
    Xw = [];
end

% negative peaks
Y = select_peaks(id,w,seq,-1,min_peak_weight,min_peak_width,peak_extend);
nn = size(Y,1);
fprintf('  negative peaks: %d peaks identified\n',nn);
if (nn > 0)
    Yu = unique_peaks(Y(:,[1 5]));
    fprintf('  negative peaks: %d unique peaks identified (%.1f%%)\n',size(Yu,1),100*size(Yu,1)./nn);
    plen = cellfun(@length,Yu(:,2));
    pscore = abs(cell2mat(Yu(:,1))./plen);
    fprintf('    mean length = %.1f (min=%d, max=%d)\n', mean(plen),min(plen),max(plen));
    fprintf('    mean absolute height = %.2f (min=%.2f, max=%.2f)\n', mean(pscore),min(pscore),max(pscore));
    Yw = unique_peak_sequences(Yu);
    fprintf('  negative peaks: %d unique peak sequences identified (%.1f%%)\n\n', size(Yw,1),100*size(Yw,1)./nn);
else
    Yu = [];
    Yw = [];
end

% combined
Z = [X;Y];
Zu = [Xu;Yu];
Zw = [Xw;Yw];
if (~isempty(Z))
    Z = sortrows(Z,-1);
    Zu = sortrows(Zu,-1);
    Zw = sortrows(Zw,-1);
end

n = size(Z,1);
if (~isempty(Z))
    fprintf('  Total: %d peaks identified (%.0f%% positive, %.0f%% negative)\n',n,100*np./n,100*nn./n);
    fprintf('  Total: %.1f%% of bases are under peaks\n', 100*sum(cellfun(@length,Z(:,5)))./numel(seq));
    fprintf('  Total: %d unique peaks identified (%.1f%%)\n',size(Zu,1),100*size(Zu,1)./n);
    fprintf('  Total: %d unique peak sequences identified (%.1f%%)\n\n',size(Zw,1),100*size(Zw,1)./n);
end

       



function X = select_peaks(id,w,seq,peak_sign,min_peak_weight,min_peak_width,peak_extend)
% find peaks in input weights

rid = repmat((1:size(w,1))',1,size(w,2))';
cid = repmat((1:size(w,2)),size(w,1),1)';

mins = find(sum(abs(w),1)>0,1,'first');
maxe = find(sum(abs(w),1)>0,1,'last');

if (peak_sign > 0)
    q = find(w' > min_peak_weight)'; % positive peaks
else
    q = find(-1*w' > min_peak_weight)'; % negative peaks
end

X = [];
if (isempty(q))
    return;
end

s = q([0 q(2:end) == q(1:end-1)+1]==0)';
e = q([q(2:end) == q(1:end-1)+1 0]==0)';
if (sum(rid(s)~=rid(e))==0)
    
    % extend peaks
    P = [rid(s) cid(s)-peak_extend cid(e)+peak_extend];
    P(P(:,2)<mins,2) = mins;
    P(P(:,3)>maxe,3) = maxe;
    
    % merge overlapping peaks
    %i = [0; (P(1:end-1,1) == P(2:end,1)).*(P(1:end-1,3) >= P(2:end,2))];
    %j = find(i==1);
    %P(j-1,3) = P(j,3);
    %P = P(i==0,:);
    
    m = size(P,1);
    S = cell(m,1);
    W = zeros(m,1);
    for k = 1:m
        S{k} = seq(P(k,1),P(k,2):P(k,3));
        W(k) = sum(w(P(k,1),P(k,2):P(k,3)));
    end
    X = [num2cell(W) id(P(:,1)) num2cell(P(:,2:end)) S];

    % filter by peak length & score
    l = cellfun(@length,S);
    if (peak_sign > 0) % positive peaks
        j = (l >= min_peak_width+2*peak_extend).*(W./l >= 0.5*min_peak_weight); 
    else % negative peaks
        j = (l >= min_peak_width+2*peak_extend).*(-1*W./l >= 0.5*min_peak_weight);
    end
    fprintf('filtering %d of %d peaks by weight & length\n', sum(j==0),size(X,1));
    X = X(j==1,:);
end

function U = unique_peaks(Z)
% U = [score] [sequence]

W(:,1) = cell2mat(Z(:,1));
[u,~,t] = unique(Z(:,2));
W(:,2) = t;
W = unique(W,'rows');

U = [num2cell(W(:,1)) u(W(:,2))];
U = sortrows(U,-1);

function U = unique_peak_sequences(Z)
% U = [score] [sequence]

[W1,~,t] = unique(Z(:,2));
W2 = accumarray(t,cell2mat(Z(:,1)));%,[],@mean);
U = [num2cell(W2) W1];
U = sortrows(U,-1);

