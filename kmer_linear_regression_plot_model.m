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

function kmer_linear_regression_plot_model(W,w0,K,n)
% W = weights (numbers)
% K = kmer ids (strings)

if (nargin < 4)
    n = 70;
end

if (~isempty(W))
    S = max(abs(W));%sum(abs(W));
    F = W./S;
    
    Q = sortrows([K num2cell(F)],[2,1]);
    q = cell2mat(Q(:,2));
    
    % select active weights
    Q = Q(abs(q)>0,:);
    q = cell2mat(Q(:,2));
    mx = max(abs(q));
    
    % select top weights
    s = sort(abs(q));
    if (max(size(s)) < n)
        n = max(size(s));
    end
    w = abs(q)>=s(end-n+1);
    n = sum(w);
    
    % plot
    barh(1:n,q(w),'FaceColor',[0.8 0.8 0.8]);
    set(gca,'ytick',1:n,'yticklabel',Q(w,1),'fontname','courier','fontsize',8);
    axis tight;
    set(gca,'xlim',[-1*mx mx]);
    xlabel('weights');
    title(sprintf('Model: Y = %.2f + %.1f * W^TX, top %d w_i out of %d\n',w0,S,n,max(size(q))));
else
    set(gca,'fontname','courier','fontsize',8);
    title(sprintf('Model: Y = %.2f\n',w0));
end
