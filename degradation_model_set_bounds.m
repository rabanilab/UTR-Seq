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

function [a,b,t] = degradation_model_set_bounds(a0,b0,t0,tm0)

if (nargin < 4)
    tm0 = 1;
end

%%% degradation rate bounds
LBa = -1*log(2)./24;     % hl > 24 hours (dg ~ 0.03)
UBa = -1*log(2)./0.25;   % hl < 15 min   (dg ~ 2.77)

a = a0;
k = a>LBa;
a(k) = LBa;
k = a<UBa;
a(k) = UBa;


%%% logX0 bounds
if (nargin > 1)
    LBb = -20;
    UBb = 20;
    
    b = b0;
    k = b<LBb;
    b(k) = LBb;
    k = b>UBb;
    b(k) = UBb;
else
    b = [];
end


%%% t0 bounds
if (nargin > 2)
    t = t0;
    t(max(a,[],2)>=LBa) = tm0;
else
    t = [];
end


% % final
% dg = -1*a;
% hl = log(2)./dg;
