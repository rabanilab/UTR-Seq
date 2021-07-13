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

function [r,df] = kmer_linear_regression_cont(data_files, feature_files, lasso_file, fit_param, result_dir)
% run kmer linear regression using given selected features
% (by lasso regularized regression)
%
% Input:
%  data_files = names of files with input data (fitted degradation model)
%  feature_files = names of files with features (e.g., kmers) for linear regression
%  lasso_file = file with lasso-regularized selected kmers
%  fit_param = parameter to fit: 0=dg, 1=t0, 2=dg with linear regression
%  ouptut_name = name of output files
%  result_dir = to save results
%  I = use to pre-select inputs (one column, for ALL input files)
%
% Output:
%  r = r-squared between input param and fitted param, per fitted model
% df = number of degrees-of-freedom in fitted model


[idx,D,tm,logx,~,bx,tx] = load_regression_data(data_files,fit_param);
[K,Krows] = load_kmer_counts(idx,feature_files,lasso_file);

%%% linear regression
file_name = [result_dir '/run_linear.out.mat'];
if (exist(file_name, 'file') ~= 2)
    
    if (fit_param == 3)
        [X,Y] = linear_fit_data(K',logx,bx,tx,tm);
    else
        X = K';
        Y = D;
    end
    [B0,b0] = linear_fit(X,Y,2);
    save(file_name,'Krows','B0','b0','-v7.3');
else
    load(file_name);
end

% model predictions
file_name = [result_dir '/run_linear.cv.mat'];
if (exist(file_name, 'file') ~= 2)
    D0 = K'*B0 + b0;
    Dj = cross_validation(D,K,logx,bx,tx,tm,fit_param==3);
    save(file_name,'idx','D0','Dj','D','-v7.3');
else
    load(file_name);
end
df = size(K,1);
r(1,1) = corr(D0,D);
r(1,2) = corr(Dj,D);
fprintf('fit: r = %.2f, rcv = %.2f (%d df)\n', r(1,1),r(1,2),df(1));

% plots
if (fit_param == 1)
    g = D>1;
    
    u = unique(D0);
    mint = 0;
    mine = inf;
    for i = 1:max(size(u))
        c = zeros(size(g));
        c(D0>=u(i)) = 1;
        e = sum(c~=g);
        if (e < mine)
            mine = e;
            mint = i;
        end
    end
    g0 = zeros(size(g));
    g0(D0>=u(mint)) = 1;
    fprintf('class: thr = %f, e = %.2f%% (%d,%d)\n',u(mint),100*sum(g0~=g)./size(g,1),sum((g0==1).*(g==0)),sum((g0==0).*(g==1)));
    plot_model_fit_t(D,D0,g,g0,size(B0,1),[result_dir '/linear']);
    
    u = unique(Dj);
    mint = 0;
    mine = inf;
    for i = 1:max(size(u))
        c = zeros(size(g));
        c(Dj>=u(i)) = 1;
        e = sum(c~=g);
        if (e < mine)
            mine = e;
            mint = i;
        end
    end
    gj = zeros(size(g));
    gj(Dj>=u(mint)) = 1;
    fprintf('class cv: thr = %f, e = %.2f%% (%d,%d)\n',u(mint),100*sum(gj~=g)./size(g,1),sum((gj==1).*(g==0)),sum((gj==0).*(g==1)));
    plot_model_fit_t(D,Dj,g,gj,size(B0,1),[result_dir '/linear_cv']);
else
    plot_model_fit_dg(D,D0,size(B0,1),[result_dir '/linear']);
    plot_model_fit_dg(D,Dj,size(B0,1),[result_dir '/linear_cv']);
    %plot_data_fit(tm,logX,b,t,D,D0,size(B0,1),[result_dir '/linear']);
    %plot_data_fit(tm,logX,b,t,D,Dj,size(B0,1),[result_dir '/linear_cv']);
end

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
kmer_linear_regression_plot_model(B0,b0,Krows,50);
saveas(h, [result_dir '/linear_model.jpg'],'jpg');

S = sum(abs(B0));
F0 = B0./S;
write_text_file([result_dir '/linear_model.txt'], ...
    sortrows([Krows num2cell([B0 F0])],[2,1]));

close all;


function [idx,D,tm,logx,ax,bx,tx] = load_regression_data(data_files,fit_param)

nfiles = max(size(data_files));
ax = [];
bx = [];
tx = [];
rx = [];
idx = [];
logx = [];
for j = 1:nfiles
    [a,b,t,~,logY,~,~,~,id,~,r,X] = degradation_model_load(data_files{j});
    tm = X{1}{1};
    ax = [ax; a];
    bx = [bx; b];
    tx = [tx; t];
    rx = [rx; r];
    idx = [idx; id];
    if (fit_param == 3)
        logx = [logx; logY];
    else
        logx = [logx; zeros(size(b))];
    end
end

if (fit_param == 1)
    D = tx;
else
    D = -1*ax(:,end);
end



function [K,Krows] = load_kmer_counts(idx,feature_files,lasso_file)

load(lasso_file,'i','Krows');
Xrows = Krows(i);

K = [];
Krows = [];
for j = 1:size(feature_files,1)
    Kj = [];
    Kjcol = [];
    for i = 1:size(feature_files,2)
        [Kij, Kjrows, Kijcol] = kmer_data_load_file([],feature_files{j,i});
        Kj = [Kj Kij];
        Kjcol = [Kjcol Kijcol];
    end
    Kjrows = regexprep(Kjrows,'_','.');
    x = ismember(Kjcol',idx);
    Kjcol = Kjcol(x);
    Kj = Kj(:,x);
    
    if ((max(size(Kjcol)) ~= max(size(idx))) || (sum(strcmp(Kjcol,idx))<max(size(idx))))
        fprintf('error loading kmer file %s\n', feature_files{j,1});
    else
        K = [K; Kj];
        Krows = [Krows; Kjrows];
    end
end

[~,j] = ismember(Krows,Xrows);
if (sum(j>0) < size(Xrows,1))
    fprintf('Error: missing input features\n');
else
    K = K(j>0,:);
    Krows = Krows(j>0);
end

function [X,Y] = linear_fit_data(K,logX,logX0,t0,tm)
% extract data for linear fit:
% log X0/Xt = t*b0 + t*sum{Kij*Bj}   t > t0
%             0                      t <=t0
% X = K
% Y = 1/t * log X0/Xt

[n,m] = size(logX);
T = repmat(tm,n,1) - repmat(t0,1,m);
T(T<0) = 0;
logR = (repmat(logX0,1,m) - logX);

X = repmat(K,m,1);
Y = logR(:)./T(:);
i = T(:)>0;

X = X(i,:);
Y = Y(i);

function [Dj,B0,b0] = cross_validation(D,K,logX,b,t,tm,fit_data_type,n)

if (nargin < 8)
    n = 10;
end

m = size(D,1);
RND(:,1) = (1:m)';
RND(:,2) = randi(n,m,1);

Dj = zeros(m,1);
b0 = zeros(n,1);
B0 = zeros(n,size(K,1));
S = zeros(n,3);
for i = 1:n
    ki = RND(RND(:,2)~=i,1);
    kj = RND(RND(:,2)==i,1);
    
    if (fit_data_type > 0)
        [X,Y] = linear_fit_data(K(:,ki)',logX(ki,:),b(ki),t(ki),tm);
    else
        X = K(:,ki)';
        Y = D(ki);
    end
    [B0(i,:),b0(i)] = linear_fit(X,Y,2);
    dj = K(:,kj)'*B0(i,:)' + b0(i);
    S(i,:) = [i size(kj,1) corr(D(kj),dj)];
    num2cell(S)
    Dj(kj) = dj;
end


function plot_model_fit_dg(D,D0,df,outpref)

dlim = [0 1];

h = figure;
hold on;
dscatter(D,D0);
line(dlim,dlim,'LineStyle', '--','color','k');
hold off;
axis tight;
axis square;
xlabel('direct estimation (rate, 1/hr)');
ylabel('regression model (rate, 1/hr)');
set(gca,'xlim',dlim,'ylim',dlim);
title(sprintf('df = %d, n = %d, r = %.2f', df,size(D0,1),corr(D0,D)));
saveas(h, [outpref '.dg.jpg'],'jpg');

close all;

function plot_model_fit_t(T,T0,M,M0,df,outpref)

tlim = [0 7];
clim = [0 1];

h = figure;
hold on;
dscatter(T,T0,'MSIZE',15);
line(tlim,tlim,'LineStyle', '--','color','k');
hold off;
axis tight;
axis square;
xlabel('direct estimation');
ylabel('lasso linear regression');
set(gca,'xlim',tlim,'ylim',tlim);
title(sprintf('df = %d, r = %.2f', df, corr(T,T0)));
saveas(h, [outpref '.t0.jpg'],'jpg');

clf;
u = unique([M(:);M0(:)]);
for i = 1:max(size(u))
    y(i,:) = hist(M(M0 == u(i)),u);
end
bar((y./sum(y(:)))','stacked');%repmat(sum(y),size(y,1),1))','stacked');
axis tight;
axis square;
set(gca,'xticklabel',strcat('class',cellstr(num2str(u))));
ylabel('fraction');
set(gca,'ylim',clim);
legend(strcat('classified as ',cellstr(num2str(u))),'location','bestOutside','box','off')
title(sprintf('df = %d, s = %.1f%%', df, 100*sum(M==M0)./size(M,1)));
box off;
saveas(h, [outpref '.c.jpg'],'jpg');





function plot_data_fit_dg(tm,logX,b,t,D,D0,df,outpref)

xlim = [0 10];
logY = degradation_model(tm,-1*D,b,t);
logY0 = degradation_model(tm,-1*D0,b,t);

h = figure;
subplot(1,2,1);
hold on;
dscatter(logX(:),logY(:));
line(xlim,xlim,'LineStyle', '--','color','k');
hold off;
axis tight;
axis square;
xlabel('measurements (log RNA)');
ylabel('model (log RNA)');
set(gca,'xlim',xlim,'ylim',xlim);
title(sprintf('r = %.2f', corr(logX(:),logY(:))));
subplot(1,2,2);
hold on;
dscatter(logX(:),logY0(:));
line(xlim,xlim,'LineStyle', '--','color','k');
hold off;
axis tight;
axis square;
xlabel('measurements (log RNA)');
ylabel('regression model (log RNA)');
set(gca,'xlim',xlim,'ylim',xlim);
title(sprintf('df = %d, r = %.2f', df, corr(logX(:),logY0(:))));
saveas(h, [outpref '.data.jpg'],'jpg');

function plot_data_fit_t(tm,logX,b,a,T,T0,df,outpref)

logY = degradation_model(tm,a,b,T);
logY0 = degradation_model(tm,a,b,T0);

h = figure;
lim = [0 10];
subplot(1,2,1);
hold on;
dscatter(logX(:),logY(:));
line(lim,lim,'LineStyle', '--','color','k');
hold off;
axis tight;
axis square;
xlabel('measurements (log RNA)');
ylabel('model (log RNA)');
set(gca,'xlim',lim,'ylim',lim);
title(sprintf('r = %.2f', corr(logX(:),logY(:))));
subplot(1,2,2);
hold on;
dscatter(logX(:),logY0(:));
line(lim,lim,'LineStyle', '--','color','k');
hold off;
axis tight;
axis square;
xlabel('measurements (log RNA)');
ylabel('regression model (log RNA)');
set(gca,'xlim',lim,'ylim',lim);
title(sprintf('df = %d, r = %.2f', df, corr(logX(:),logY0(:))));
saveas(h, [outpref '.data.jpg'],'jpg');


