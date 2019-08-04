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

function [score,opt_pos,opt_seq] = peaks_pwm_score(seq,PWM,pwmC,pseudoC)
% calculate log-ratio score by given PWM matrix
%   score = P(seq|PWM)

if (nargin < 3)
    pwmC = 20;
end
if (nargin < 4)
    pseudoC = 1; % 1/25 ~ 4%
end



% reassign PWMs with pseudo-counts
m = size(PWM,2);
log2pwm = add_pseudocounts(PWM,pwmC,pseudoC);

% score sequences
k = size(seq,1);
score = zeros(k,1);
opt_pos = zeros(k,1);
opt_seq = cell(k,1);
for i = 1:k
    seqi = [repmat('-',1,m) seq{i} repmat('-',1,m)];
    max_score = -inf;
    max_pos = 0;
    for istart = 1:(length(seqi)-m)
        iend = istart + m - 1;
        scorei = get_score(seqi(istart:iend),log2pwm);
        if (scorei > max_score)
            max_score = scorei;
            max_pos = istart;
        end
    end
    score(i) = 2.^(max_score);
    opt_pos(i) = max_pos - m;
    opt_seq{i} = seqi(max_pos:(max_pos+m-1));
end


function score = get_score(seq,log2pwm)

L = 'ACGT-';
score = 0;
m = size(log2pwm,2);
for j = 1:m
    score = score + log2pwm(L==seq(j),j);
end

function log2pwm = add_pseudocounts(PWM,pwmC,pseudoC)

[n,m] = size(PWM);
C = round(pwmC.*PWM);
C = C + pseudoC;
C = [C; ones(1,m)];
log2pwm = log2(C) - repmat(log2(sum(C)),n+1,1);
