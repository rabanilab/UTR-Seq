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

function [W,Wid,xPWM,xPWMstat,xPWMw] = peaks_assign_seeds(pseq,pscore,bgseq,PWM,PWMstat,assign_type)
% Output:
%  W = list of peaks [peak seq] [peak score] [motif score] [pos/neg] [assigned motif]
%  PWMid = motif ids
%  PWM = motif PWMs
%  PWMstat = peak statistics

if (nargin < 6)
    assign_type = 1;
end


% background distribution
[a,~,t] = unique([bgseq{:}]);
c = accumarray(t,1);
BG = c./sum(c);
%[cellstr(a') num2cell(BG)]

% assign peaks
PWMid = PWMstat(:,1);
PWMn = cell2mat(PWMstat(:,2));
PWMs = cell2mat(PWMstat(:,3));

[W1,p1] = assign_peaks(BG,pseq(pscore>0),pscore(pscore>0),PWMs,PWMid,PWM,1);
%fprintf('Assign %d positive PWM models\n', sum(p1));
%[Q(p1,[1 3 2]) num2cell(100*PWMn(p1)./sum(PWMn))]

[W2,p2] = assign_peaks(BG,pseq(pscore<0),pscore(pscore<0),PWMs,PWMid,PWM,-1);
%fprintf('Assign %d negative PWM models\n', sum(p2));
%[Q(p2,[1 3 2]) num2cell(100*PWMn(p2)./sum(PWMn))]

W = [W1;W2];
p = p1+p2 > 0;

Wid = strcat(W(:,end),':',num2str(cell2mat(W(:,end-1))+1));
Wid = regexprep(Wid,':2',':P');
Wid = regexprep(Wid,':0',':N');
[u,~,t] = unique(Wid);
c = accumarray(t,1);
w = accumarray(t,abs(cell2mat(W(:,2))));
[u1,u2] = strtok(u,':');
xPWMstat = [u1 num2cell([c 2*cellfun(@isempty,regexp(u2,':N'))-1])];
xPWMw = w;

xPWM = PWM(p);
xPWMid = strcat(PWMid(p),':',num2str(PWMs(p)+1));
xPWMid = regexprep(xPWMid,':2',':P');
xPWMid = regexprep(xPWMid,':0',':N');

[i,j] = ismember(xPWMid,u);
if (assign_type == 0)
    ts = xPWMstat;
    tw = xPWMw;
    [u1,u2] = strtok(xPWMid,':');
    xPWMstat = [u1 num2cell([zeros(size(xPWMid,1),1) 2*cellfun(@isempty,regexp(u2,':N'))-1])];
    xPWMw = zeros(size(xPWMid,1),1);
    xPWMstat(i,:) = ts(j(j>0),:);
    xPWMw(i,:) = tw(j(j>0),:);
else
    xPWM = xPWM(i);
    xPWMstat = xPWMstat(j(j>0),:);
    xPWMw = xPWMw(j(j>0));
end



function [W,p] = assign_peaks(BG,pseq,pscore,PWMs,PWMid,PWM,d)

n = size(pseq,1);
if (d>0)
    p = PWMs>0;
else
    p = PWMs<0;
end

W = [];
if (n>0)
    [S,N] = pwm_select(pseq,PWM(p),PWMid(p),BG);
    S(N<0) = {'NA'};
    W = [W;[pseq num2cell([pscore N d*ones(n,1)]) S]];
end

function [ids,score,pos] = pwm_select(seq,PWM,pwm_ids,BG)
% assign peaks into best fitting PWM
%   score = log2(P(seq|PWM)) - log2(P(seq|BG))

if (nargin < 4)
    BG = [0.25; 0.25; 0.25; 0.25];
end

n = max(size(PWM));
S = ones(size(seq,1),n);
for i = 1:n
    l = size(PWM{i},2);
    bg = log2(peaks_pwm_score(seq,repmat(BG,1,l)));
    S(:,i) = log2(peaks_pwm_score(seq,PWM{i}));
    S(:,i) = S(:,i) - bg;
end

[score,pos] = max(S,[],2);
ids = pwm_ids(pos);
