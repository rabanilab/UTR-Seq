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

function utrseq_fit_deg_rates(tm, data_file_name, result_dir)
% logX = log2 normalized counts per time
%    C = raw counts per time
% output_name = '3U.00A.seq1022';

utm  = 1:0.1:10;
alpha = 0.05;
error_func = 0;

if (exist(result_dir,'dir') ~= 7)
    mkdir(result_dir);
end

%%% load data
fprintf('** INPUT DATA\n');
[cid,Nid,N,Qid,Q] = utrseq_load_counts(data_file_name);
fprintf('load %s (%d x %d)\n', data_file_name, size(Nid,1), size(cid,2));
fprintf('\n');

%%% normalize
si_cnt = 50;
si_frac = 0.001;
si_type = 0;

fprintf('** Normalize DATA\n');
[~,~,N,logX,Qr,Qm,Qn] = utrseq_normalize({1:size(N,2)},Nid,cid,tm,Q,N,si_cnt,si_frac,si_type);
logX = logX{1};
N = N{1};
save([result_dir '/degradation.norm.mat'],'cid','tm','Qid','Q','Nid','N','logX','Qr','Qm','Qn');
write_text_file([result_dir '/degradation.norm.txt'],[Nid num2cell(logX)]);

[I,O,logZ,iC] = samples_select({tm},{N},{logX},utm);
O = O{1};
logZ = logZ{1};
logA = samples_min(iC,{logZ},5);

S = samples_variation({logZ});
S = S{1};
fprintf('\n');

%%% optimization
fprintf('** Optimization\n');
[a0,b0,t0,r0] = degradation_model_fit(tm(O),logZ,logA,S,error_func);
fprintf('r-squared test:\n');
fprintf(' r-squared >= 0.7, n=%d (%.1f%%)\n', sum(r0>=0.7),100*sum(r0>=0.7)./size(r0,1));
fprintf(' r-squared >= 0.9, n=%d (%.1f%%)\n', sum(r0>=0.9),100*sum(r0>=0.9)./size(r0,1));

p0 = degradation_model_goodnessoffit(tm(O),logZ,S,logA,a0,b0,t0);
fprintf('goodness of fit test:\n');
fdr_alpha = fdr(alpha,p0);
fprintf(' reject model at %.0f%% FDR, n=%d (%.1f%%)\n', 100*alpha, sum(p0<=fdr_alpha),100*sum(p0<=fdr_alpha)./size(p0,1));
fprintf('\n');

%%% save data to output files
m0 = size(Nid,1);
A = nan(m0,size(a0,2));
A(I,:) = a0;
B = nan(m0,1);
B(I) = b0;
T = nan(m0,1);
T(I) = t0;
R = nan(m0,1);
R(I) = r0;
X = nan(m0,1);
X(I) = 1;
save_data([result_dir '/degradation.fit'],Nid,tm(O),A,B,T,R,X,S,logX(:,O),N(:,O),logA);


function [Cid,Sid,S,Qid,Q] = utrseq_load_counts(data_file_name)

X = importdata(data_file_name);
% f = fopen(data_file_name);
% Z = textscan(f,'%s %s %s %s %s %s %s %s %s %s %s %s');
% fclose(f);
% X.textdata = [Z{:}];
% X.data = str2double(X.textdata(2:end,2:end));

% counts
S = X.data;
Sid = X.textdata(2:end,1);
Cid = X.textdata(1,2:end);

% spike-ins counts
k = strncmp(Sid,'SPIKE_',6);
Q = S(k,:);
Qid = Sid(k);
S = S(k==0,:);
Sid = Sid(k==0,:);


function save_data(fname,id,tm,a,b,t,r,x,sd,logY,C,logA)

save([fname '.mat'],'id','tm','a','b','t','r','x','sd','logY','C','logA');
write_text_file([fname '.txt'],[id num2cell([-1*log(2)./a a b t r])]);


function [a,b,t,r] = degradation_model_fit(tm,logX,logA,S,error_func)
% Model: 
%     dX/dt = -aX
%      X(t) = X0 * exp(-at)
%  log X(t) = log X0 - at
% log2 X(t)./log2(exp(1)) = log2 X0./log2(exp(1)) - at
% log2 X(t) = log2 X0 - a*log2(exp(1)) t

% error function
if (error_func == 1)
    %fprintf('least square error\n');
    e_func = @lsq_error;
else
    %fprintf('r-square error\n');
    e_func = @rsq_error;
end

[a,b,t] = fit_degradation_2p(min(tm), tm, logX, logA, S, e_func);
[a,b,t] = degradation_model_set_bounds(a,b,t,min(tm));
logY = degradation_model(tm,a,b,t,logA);
b = b - nanmean(logY-logX,2);
[a,b,t] = degradation_model_set_bounds(a,b,t,min(tm));
r = r_func(tm,logX,logA,a,b,t);


function [a,b,t,e] = fit_degradation_2p(t0, tm, logX, logA, S, e_func)
% fit 2 gene specific parameter (b,a)

% model parameters
n = size(logX,1);
t = t0*ones(n,1);

% linear model
w = repmat(1./S,n,1);
%logX(logX<logA) = logA;
w(logX <= logA) = 0;
w(isnan(logX)) = 0;
x = repmat(tm - t0,n,1);
x(x<0) = 0;
[a,b] = linear_fit_2p(x,logX,w);

% error
e = e_func(tm,logX,logA,a,b,t,S);
%fprintf('%.1f: e = %.2e\n',t0, sum(e));


function r = r_func(tm,logX,logA,a,b,t)
% r squared

logY = degradation_model(tm,a,b,t,logA);
r = r_squared(2.^logX,2.^logY);

function e = rsq_error(tm,logX,logA,a,b,t,S)
% r squared error

logY = degradation_model(tm,a,b,t,logA);
e = 1-r_squared(2.^logX,2.^logY);

function e = lsq_error(tm,logX,logA,a,b,t,S)
% least squares error

logY = degradation_model(tm,a,b,t,logA);
e = least_squares(2.^logX,2.^logY,S);



function minX = samples_min(C,logX,minR)
% set a minimal expression value
%
% Input:
%   C = counts (ids x times)
%   logX = normalized expression
%   minR = minimal read count
% Output:
%   minX = minimal expression value for logX


if ((nargin < 3) || (isempty(minR)))
    minR = 5;
end

k = max(size(logX));
min_values = [];
for i = 1:k
    minC = min(C{i}(C{i} > minR));
    min_values = [min_values;unique(logX{i}(C{i} == minC))];
end

minX = prctile(min_values,25);%mean(min_values);
for i = 1:k
    n = sum(logX{i}(:)<minX);
    fprintf('minimal expression: minX = %.2f, %.1f%% below min (%d)\n', ...
        minX,100*n./max(size(logX{i}(:))),n);
    str = sprintf('%d, ',sum(logX{i}<minX));
    fprintf('reporters below min: (%s)\n', str(1:end-2));
    
    y = hist(sum(logX{i}<minX,2),0:size(logX{i},2));
    num2cell([(0:size(logX{i},2))' y'])
end

function [I,O,ilogX,iC,utm] = samples_select(tm,C,logX,utm,minR,minF)
% select reporters by minimal expression over temporal samples
%
% Input:
%   C = counts (ids x times)
%   logX = normalized expression
%   minR = minimal sample-average read count per instance (row)
%   minF = minimal fraction of reads per sample (column)
% Output:
%   I = 0/1 selection
%   ilogX = normalized expression of for selected samples
%   iC = counts of for selected samples

k = max(size(logX));
n = size(logX{1},1);

if ((nargin < 5) || (isempty(minR)))
    minR = 20;
end
if ((nargin < 6) || (isempty(minF)))
    minF = 0.035;
end

% read coverage per sample
mnu = min([1 utm]);
mxu = max([utm 10]);
O = cell(1,k);
for i = 1:k
    x = sum(C{i},1)./sum(C{i}(:));
    if (max(size(x))>1)
        m = ([x(1) x(1:end-2) 0]+[x(2) x(2:end-1) x(end-1)]+[0 x(3:end) x(end)])/3;
    else
        m = x;
    end
    O{i} = (tm{i}>=mnu).*(tm{i}<=mxu).*(x>=minF).*(x>=0.5*m) == 1;
    num2cell([tm{i}' x' m' O{i}'])
    s = sort(tm{i}(O{i}));
    if (max(size(s))>2)
        utm = utm(utm<=s(end-2)); % must use at least 3 points to fit a linear line
    else
        utm = utm(1);
    end
end
fprintf('sample selection threshold: frac. reads >= %.1f%%\n', 100*minF);
fprintf('samples selected: %d/%d (%.1f%%)\n', sum([O{:}]), max(size([O{:}])), 100*sum([O{:}])./max(size([O{:}])));

% read coverage per instance
M = zeros(n,k);
for i = 1:k
    if (isempty(C))
        M(:,i) = minR;
    else
        M(:,i) = sum(C{i}(:,O{i}),2)./sum(O{i},2);
    end
end
I = (sum(M >= minR,2) == k);
fprintf('reporter selection threshold: avg. reads >= %d\n', minR);
fprintf('reporters selected: %d (%.1f%%)\n', sum(I), 100*sum(I)./n);

% final
ilogX = cell(1,k);
iC = cell(1,k);
for i = 1:k
    ilogX{i} = logX{i}(I,O{i});
    if (~isempty(C))
        iC{i} = C{i}(I,O{i});
    end
end

function S = samples_variation(logX,type,F)
% calculate expression variability per sample
% mean = sum(x_i)/n
%  var = sum((x_i - mean)^2)/(n-1)
%  std = sqrt(var)
%
% Input:
%  logX = cell array of replicated samples
%         each sample is a matrix of [reporters]x[times]
%  type = how to estimate variability: 
%         0=std [default], 1=var, 2=coefficient of variation
% Output:
%     S = variability per sample (one number per time)

if ((nargin < 2) || isempty(type))
    type = 0;
end
if (nargin < 3)
    F = 0.32;
end

n = max(size(logX));
S = cell(1,n);
for i = 1:n
    x = logX{i};
    k = size(x,2);
    
    v = zeros(1,k);
    s = zeros(1,k);
    m = zeros(1,k);
    c = zeros(1,k);
    for j = 1:k
        %lb = prctile(x(:,j),10);
        %ub = prctile(x(:,j),90);
        %o = (x(:,j)>lb).*(x(:,j)<ub)==1;
        pt = 100 - ceil(1e6./size(x(:,j),1));
        lb = prctile(x(:,j),max(pt,0));
        o = (x(:,j)>lb)==1;
        v(j) = var(x(o,j),0,1,'omitnan');
        s(j) = std(x(o,j),0,1,'omitnan');
        m(j) = prctile(x(o,j),50);
        c(j) = sum(o);
    end
    m(m<=0) = min(x(x>0));
    
    if (type == 1)
        S{i} = mean(v,1);
    elseif (type == 2)
        S{i} = mean(s./m,1);
    else
        S{i} = mean(s,1);
    end
       
    str = sprintf('%.2f, ',S{i});
    fprintf('standard deviation: (%s)\n', str(1:end-2));
    mn = min(S{i});
    mx = max(S{i});
    S{i} = F*0.5*(2 + (S{i}-mn)./(mx-mn));
    str = sprintf('%.2f, ',S{i});
    fprintf('corrected standard deviation: (%s)\n', str(1:end-2));
    str = sprintf('%d, ',c);
    fprintf('number of reporters used: (%s)\n', str(1:end-2));
end


function [a,b] = linear_fit_2p(x,y,w)
% linear regression: y = a*x + b
% w = weights for fit (higher weight = stronger influence on result)

minW = 1e-10;
maxW = 1;

if (nargin < 3)
    w = ones(size(x));
end

% weights
w(abs(w)<minW) = minW;
w(abs(w)>maxW) = maxW;
w2 = w.^2;%1./(w.^2);

% linear fit
sx = nansum(x.*w2,2);
mx = nansum(x.*w2,2)./sum(w2,2);
my = nansum(y.*w2,2)./sum(w2,2);
sxx = nansum(x.*w2.*x,2);
sxy = nansum(x.*w2.*y,2);

a = (sxy - my.*sx)./(sxx - mx.*sx);
a(isnan(a)) = 0;
b = my - a.*mx;

i = sum(w>minW,2)<2;
a(i) = 0;
b(i) = nanmean(y(i,:),2);

% finalize fit
a = a./log2(exp(1)); % convert natural to log2 base
[a,b] = degradation_model_set_bounds(a,b); % boundary conditions

function r = r_squared(X,Y)
% X = measurements
% Y = model predictions
% r = r squared, coefficient of determination

[n,m] = size(X);

mn = repmat(nanmean(X,2),1,m);
s1 = nansum((X - Y).^2,2); % SSerr
s2 = nansum((X - mn).^2,2); % SStot
r = zeros(n,1);
r(s2>0) = 1 - s1(s2>0)./s2(s2>0); % assume: s1<s2
r(r<0) = 0;

function fdr_alpha = fdr(alpha, p)
% fdr

vp = p(:);
if (sum(vp>0)>0)
    vp(vp == 0) = min(vp(vp > 0));
end
m = size(vp,1);
pk = sort(vp, 'ascend');
fdr_alpha = ((1:m)')*(alpha/m);

pass_test = find(pk <= fdr_alpha);
if (isempty(pass_test))
    fdr_alpha = 0;
else
    k = max(pass_test);
    fdr_alpha = alpha*k/m;
end

function [ID,TM,Cr,Cn,Qr,Qm,Qn] = utrseq_normalize(S,rid,cid,tm,Q,C,si_min_count,si_fraction,si_norm_type)
% output data is per series:
%   ID = cids
%   TM = times
%   Qr = raw spike-in counts
%   Cr = raw UTR-Seq counts
%   Qn = normalized spike-in counts
%   Cn = normalized UTR-Seq counts

n = max(size(S));
if (nargin < 8)
    si_min_count = 50;
end
if (nargin < 9)
    si_fraction = 0.001;
end
if (nargin < 10)
    si_norm_type = 0;
end


ID = cell(n,1);
TM = cell(n,1);
Qr = cell(n,1);
Qm = cell(n,1);
Qn = cell(n,1);
Cr = cell(n,1);
Cn = cell(n,1);
fprintf('normalizing data: %d series, %d reporters\n', size(rid,1), n);

for i = 1:n
    fprintf('series %d\n', i);
    
    si = S{i};
    ID{i} = cid(si);
    TM{i} = tm(si);
    Qr{i} = Q(:,si);
    Cr{i} = C(:,si);
    Qm{i} = utrseq_normalize_spikes(tm(si),Qr{i},Cr{i},si_min_count,si_fraction,[]);
    [Cn{i},Qn{i}] = utrseq_normalize_counts(Cr{i},Qm{i},[],si_norm_type);
    close all;
end


function [Q3,Q2,Q1] = utrseq_normalize_spikes(tm,Q,C,N,M,plot_prefix)
% normalizing the first time point spike-ins by a linear fit
% removing outliar spikes
%
% Input:
% tm = times
%  Q = control gene counts
%  C = reporter counts
%  N = minimal sum of spike-in reads to use for regression fit
%  M = fraction of spike-ins out of total reads
%  plot_prefix = for QC plots during normalization
%
% Output: adjusted spike-in counts

if (nargin < 4)
    N = 1e4;%[];
end
if (nargin < 5)
    M = 0.001;%[];
end
if (nargin < 6)
    plot_prefix = [];
end

R = [1 1/2 1/4 1/8 1/16]';
T1 = 80;
T2 = 50;

% spike-ins: remove outliars
Q1 = Q;
if (size(Q1,1) == 5)
    j1 = Q1(2:5,:)./Q1(1:4,:) <= 0.25;
    j1 = [zeros(1,size(Q1,2)); j1];
    j2 = Q1(2:5,:)./Q1(1:4,:) >= 1;
    j2 = [j2; zeros(1,size(Q1,2))];
    Q1(j1.*j2==1) = 0;
    fprintf('spikes: %d outliars\n', sum(Q1(:)==0));
end

% spike-ins: adjust fraction of spike-ins out of total counts
Q2 = Q1;
j = (sum(Q2)>T2).*(tm>0).*(tm<24) == 1;
if (~isempty(N))
    j = j.*(sum(C)>N) == 1;
end
Q2(:,j==0) = R*sum(Q2(:,j==0))./sum(R);

if (~isempty(M))
    W = C(sum(C(:,j),2) >= prctile(sum(C(:,j),2),T1),j);
    V2 = Q2(:,j);
    P = sum(W(:))./sum(V2(:)) * (M/(1-M));
    Q2 = round(Q2*P);
    
    V3 = Q2(:,j);
    B1 = sum(V2(:))/sum([W(:);V2(:)]);
    B2 = sum(V3(:))/sum([W(:);V3(:)]);
    fprintf('spikes: adjust fraction from %.2f %% to %.2f %%\n', 100*B1, 100*B2);
end

% spike-ins: temporal logistic fit
Q3 = Q2;
if (~isempty(N))
    x1 = 1e5*sum(Q3)./sum([C;Q3]);
    y = log2(x1(j));
    y = [y(1) y(1) y(1) y]; % 4-fold weight to first sample
    x = tm(j);
    x = [x(1) x(1) x(1) x];
    if (sum(j) >= 5)
        s1 = logistic_fit(y,x,100,[0 0 0 0],[0.6 inf inf 24]); % 4 param
        z1 = 2.^logistic_eval(s1,tm);
    else
        s1 = polyfit(x,y,1); % 2 param
        z1 = 2.^polyval(s1,tm);
    end
    Q3 = round(Q3.*repmat(z1./x1,size(Q3,1),1));
    Q3(isnan(Q3)) = 0;
    Q3(isinf(Q3)) = 0;
end

p = sum(Q3)./sum(C);
p(p<1) = 1;
Q3 = round(Q3./repmat(p,size(Q3,1),1));
fprintf('spikes: %d fraction > 50%%\n', sum(p>1));


if (~isempty(plot_prefix))
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    
    hold on;
    x1 = 1e5*sum(Q1)./sum([C;Q1]);
    plot(tm,log2(x1),'-k','linewidth',1,'marker','.','markersize',30);
    x2 = 1e5*sum(Q2)./sum([C;Q2]);
    plot(tm,log2(x2),'-m','linewidth',1,'marker','.','markersize',30);
    x3 = 1e5*sum(Q3)./sum([C;Q3]);
    plot(tm,log2(x3),'-r','linewidth',1,'marker','.','markersize',30);
    L = [];
    if (~isempty(N))
        tx = min(tm):0.1:max(tm);
        if (sum(j) >= 5)
            y1 = logistic_eval(s1,tx);
            L = sprintf('sigmoid: (b,h0,h1,t1) = (%.1f,%.1f,%.1f,%.1f)',s1);
        else
            y1 = polyval(s1,tx);
            L = sprintf('linear: y = %.1fx + %.1f',s1);
        end
        plot(tx,y1,'-b','linewidth',1.5);
    end
    plot(tm(j==0),log2(x3(j==0)),'.b','markersize',30);
    hold off;
    axis tight;
    axis square;
    title({'spike-in normalization' L});
    legend({'original' 'normalized' 'fitted model'},'location','bestoutside','box','off');
    xlabel('time');
    ylabel('sum(spike-in reads)/sum(reads)');
    set(gca,'ylim',[0 15],'xlim',[0 24]);
    
    saveas(h, [plot_prefix '.normspikes.jpg'],'jpg');
end

% % spike-ins: linear fit of first time point spikes
%F = 75;
%I = 2;
% if (~isempty(I))
%     x = [];
%     y = [];
%     for i = 1:max(size(I))
%         j = (C(:,1)>=prctile(C(:,1),F)).*(C(:,I(i))>=prctile(C(:,I(i)),F))==1;
%         x = [x;log2(C(j,I(i)))];
%         y = [y;log2(C(j,1))];
%     end
%     k = (~isinf(x)).*(~isinf(y)).*(~isnan(x)).*(~isnan(y)) == 1;
%     %P = polyfit(x(k),y(k),1);
%     %fprintf('y = %.2fx + %.2f \n',P);
%     %P = robustfit(x(k),y(k));
%     %fprintf('y = %.2fx + %.2f \n',P(2),P(1));
%     P = glmfit(x(k),2.^y(k),'poisson','link','log');
%     fprintf('log(y) = %.2fx + %.2f \n',P(2),P(1));
%     
%     if (~isempty(plot_prefix))
%         subplot(1,2,2);
%         hold on;
%         if (sum(k) > 5)
%             dscatter(x(k),y(k));
%         end
%         mn = min([x(k);y(k)]);
%         mx = max([x(k);y(k)]);
%         %z = polyval(P,[mn mx]);
%         %plot([mn mx],z,'-k');
%         %z = P(1)+P(2)*[mn mx];
%         %plot([mn mx],z,'-k');
%         z = glmval(P,[mn mx],'log');
%         plot([mn mx],log2(z),'-k');
%         line([mn mx],[mn mx],'LineStyle', '--', 'color', 'r');
%         set(gca,'xlim',[mn mx],'ylim',[mn mx]);
%         hold off;
%         axis square;
%         xlabel('1 hpf');
%         ylabel('0 hpf');
%         title(sprintf('y = %.2fx + %.2f (n = %d)',P, sum(k)));
%     end
%     
%     J = zeros(size(Q(:,I)));
%     for i = 1:max(size(I))
%         J(:,i) = glmval(P,log2(Q(:,I(i))),'log');
%     end
%     Q(:,1) = round(mean(J,2));
%     fprintf('spikes: fit first time point spikes by y = %.2fx + %.2f (n = %d)\n',P,sum(k));
% end

function [logN,logM,Q,P,E] = utrseq_normalize_counts(C,Q,plot_prefix,norm_type)
% estimate a linear regression model from input spike-inns, and
% use that to normalize input counts
%
% Input:
% C = reporter counts
% Q = spike-in counts
% plot_prefix = for QC plots
% norm_type = 0/1/2
%   0: linear regression fit of spike-ins with a binomial model (default)
%   1: linear regression fit of spike-ins with a poisson model
%   2: normalized sum of spike-ins
%
% Output is log2 transformed:
% logN = normalized C
% logM = normalized Q
% P = polynomials, P = [a,b] : y = ax+b, one row per X column
% E = polynomials fit error (lsq)

if (nargin < 3)
    plot_prefix = [];
end
if (nargin < 4)
    norm_type = 0;
end

[n,m] = size(C);
fprintf('input: %d x %d matrix\n', n, m);


% calculate normalization weights
Y = [50 25 12.5 6.25 3.125 0]';
Yweights = [4 4 2 2 1 10];%[4 4 2 2 1 1]
if (norm_type == 2)
    Y = ones(size(Q,1)+1,1);
    Yweights = Y;
end

P = zeros(m,2); % [b,a] : y = ax+b
for i = 1:m
    X = [Q(:,i); 0];
    j = (X==0).*(Y>0)==0;
    if (norm_type == 1)
        P(i,:) = glmfit(Y(j),X(j),'poisson','link','identity','weights',Yweights(j));
    elseif (norm_type == 2)
        P(i,1) = 0;
        P(i,2) = nanmean(Q(:,i),1)./nanmean(Q(:));
    else
        K = sum([Q(:,i);C(:,i)]);
        P(i,:) = K*glmfit(Y(j),[X(j) K*ones(size(X(j)))],'binomial','link','identity','weights',Yweights(j));
    end
end
P(P(:,2) == 0,2) = 1;
fprintf('linear regression parameters:\n');
[{'sample' 'a' 'b'};num2cell([(1:m)' 1./P(:,2) -1*P(:,1)./P(:,2)])]

% normalize read counts
N = zeros(n,m);
M = zeros(size(Q,1)+1,m);
E = zeros(m,2);
mn0 = zeros(m,1);
mn1 = zeros(m,1);
for i = 1:m
    X = C(:,i); % counts
    N(:,i) = (X - P(i,1))./P(i,2);
    mn0(i) = min(N(C(:,i)>0,i),[],1);
    mn1(i) = min(N(N(:,i)>0,i),[],1);
    N(C(:,i)==0,i) = mn0(i)/2;
    N(N(:,i)<=0,i) = mn1(i)/2;
    
    X = [Q(:,i); 0]; % spike-ins
    M(:,i) = (X - P(i,1))./P(i,2);
    E(i,1) = sum((M(:,i)-Y).^2);
    M(X==0,i) = mn0(i)/2;
    M(M(:,i)<=0,i) = mn1(i)/2;
    E(i,2) = sum((M(:,i)-Y).^2);
end
fprintf('normalization:\n');
[{'sample' 'min0' 'min1' 'minN' 'Err0' 'Err1'};num2cell([(1:m)' mn0 mn1 min(N,[],1)' E])]

% logscale conversion
logN = log2(N);
logM = log2(M(1:end-1,:));

% plots
if (~isempty(plot_prefix))
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    
    clf;
    logC = log2(C);
    logC(isnan(logC)) = 0;
    logC(isinf(logC)) = 0;
    mn = min([min(logC(:)) min(logN(:))]);
    mn = ceil(mn);
    mx = max([max(logC(:)) max(logN(:))]);
    mx = floor(0.9*mx);
    x = mn:0.25:mx;
    for i = 1:min(m,15)
        subplot(3,5,i);
        hold on;
        y = hist(logC(:,i),x);
        plot(x,y,'-k','linewidth',1.2); % unnormalized = black
        y = hist(logN(:,i),x);
        plot(x,y,'-r','linewidth',1.2); % normalized = red
        xlabel('log2 expression');
        ylabel('number of fragments');
        axis tight;
        axis square;
        set(gca,'ylim',[0 12000],'ytick',[],'xtick',[]);
        if (i==1)
            title({sprintf('%d', i); sprintf('black=orig,red=norm')});
        else
            title(sprintf('%d', i));
        end
    end
    saveas(h, [plot_prefix '.hist.jpg'],'jpg');
    
    clf;
    xlim = [0 max(Q(:))];
    ylim = [0 50];
    for i = 1:m
        subplot(4,4,i);
        hold on;
        a = min(Q(:,i)):max(Q(:,i));
        b = (a - P(i,1))./P(i,2);
        plot(a,b,':','color','k','linewidth',1.2);
        a = prctile(C(:,i),5):prctile(C(:,i),99.99);
        b = (a - P(i,1))./P(i,2);
        plot(a,b,'-','color','k','linewidth',1.2);
        %dscatter(C(:,i),N(:,i),'MSIZE',50);
        plot([Q(:,i); 0],Y,'.r','markerSize',20);
        xlabel('reads (count)');
        ylabel('RNA (femtograms)');
        axis tight;
        axis square;
        set(gca,'xtick',xlim(1):1000:xlim(2),'ytick',sort(unique(Y)));%,'xlim',xlim,'ylim',ylim);
        title(sprintf('y=%.1fx%+.1f (e=%.1f)',1./P(i,2),-1*P(i,1)./P(i,2),E(i,2)));
    end
    saveas(h, [plot_prefix '.norm.jpg'],'jpg');
    
    clf;
    xlim = [0 15];
    ylim = [-7 11];
    for i = 1:m
        subplot(4,4,i);
        hold on;
        %x1 = log2(C(:,i));
        %x1(C(:,i)<1) = 0;
        %dscatter(x1,logN(:,i),'MSIZE',50);
        y1 = log2(Q(:,i));
        y1(Q(:,i)<1) = 0;
        y2 = log2(Y(1:end-1));
        plot(y1,y2,'.r','markerSize',25);
        a = ceil(max(P(i,1),1)):10:5000;
        b = (a - P(i,1))./P(i,2);
        plot(log2(a),log2(b),'-','color','k','linewidth',1.2);
        hold off;
        xlabel('reads (log count)');
        ylabel('RNA (log femtograms)');
        axis tight;
        axis square;
        set(gca,'xtick',-100:4:100,'ytick',unique(round(y2)),'xlim',xlim,'ylim',ylim);
        title(sprintf('y=%.1fx%+.1f (e=%.1f)',1./P(i,2),-1*P(i,1)./P(i,2),E(i,2)));
    end
    saveas(h, [plot_prefix '.lognorm.jpg'],'jpg');
end

