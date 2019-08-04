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

function [a,b,t,x,logY,logD,logH,I,id,p,r,X] = degradation_model_load(input_name,ids)
% Load model fit parameters
% [a,b,t,x,logY,logD,logH,I,id,p,r,X] = degradation_model_load(name,ids)
%
% X = {{tm logY a b t r}};

if (nargin < 2)
    ids = [];
end

lbd = -5;  % ~ log2(log(2)./24)
ubd = 1.5; % ~ log2(log(2)./0.25)

load(input_name,'a','b','t','x','r','tm','sd','logY','logA','id');
if (~isempty(ids))
    [~,i] = ismember(ids,id);
    id = id(i(i>0));
    
    ai(i>0,:) = a(i(i>0),:);
    ai(i==0,:) = NaN;
    a = ai;
    bi(i>0,:) = b(i(i>0),:);
    bi(i==0,:) = NaN;
    b = bi;
    ti(i>0,:) = t(i(i>0),:);
    ti(i==0,:) = NaN;
    t = ti;
    xi(i>0,:) = x(i(i>0),:);
    xi(i==0,:) = NaN;
    x = xi;
    ri(i>0,:) = r(i(i>0),:);
    ri(i==0,:) = NaN;
    r = ri;
    if (iscell(logY))
        n = max(size(logY));
        for j = 1:n
            logYj = logY{j}(i(i>0),:);
            logYj(i==0,:) = NaN;
            logY{j} = logYj;
        end
    else
        logYi = logY(i(i>0),:);
        logYi(i==0,:) = NaN;
        logY = logYi;
    end
end
fprintf('Load %s: %d samples\n',input_name,sum(~isnan(x),1));

t = round(10*t)/10;
if (iscell(logY))
    p = degradation_model_goodnessoffit([tm{:}],[logY{:}],[sd{:}],logA,a,b,t);
    p(isnan(b)) = nan;
else
    p = degradation_model_goodnessoffit(tm,logY,sd,logA,a,b,t);
    p(isnan(b)) = nan;
end
if (size(a,2) < 2)
    a = [zeros(size(a)) a];
    a(sum(isnan(a),2)>0,:) = nan;
end

% logD = [logdg(1h), logdg(10h)]
i = t>1;
logD(i==1,1) = log2(-1*a(i==1,1)); % t>1
logD(i==0,1) = log2(-1*a(i==0,2)); % t<=1
i = t>10;
logD(i==1,2) = log2(-1*a(i==1,1)); % t>10
logD(i==0,2) = log2(-1*a(i==0,2)); % t<=10

logD(isinf(logD)) = lbd;
logD(imag(logD)~=0) = lbd;
logD(logD<lbd) = lbd;
logD(logD>ubd) = ubd;

% H = log(2)/D
logH = log2(log(2)) - logD;

% select 'good fit' examples
if (iscell(logY))
    n = max(size(logY));
    I = (max(b,[],2) > 1.0-n*0.1).*((p > 0.05)+(r > 0.7-n*0.05)>0);
else
    I = (b > 1.0).*((p > 0.05)+(r > 0.7)>0);
end
fprintf('Load %s: %d selected samples (%.1f%%)\n',input_name,sum(I,1),100*sum(I,1)/sum(~isnan(x),1));

% a collection of all data (for plots)
if (iscell(logY))
    n = max(size(logY));
    X = cell(1,n);
    for j = 1:n
        X{j} = {tm{j} logY{j} a b t r};
    end
else
    X = {{tm logY a b t r}};
end
