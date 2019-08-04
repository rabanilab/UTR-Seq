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

function [r,df] = kmer_lasso_regression_cont(data_files, feature_files, fit_param, result_dir)
% run kmer linear regression using LASSO regularization to select
% k-mer features for prediction
%
% Input:
%  data_files = names of files with input data (fitted degradation model)
%  feature_files = names of files with features (e.g., kmers) for linear regression
%  fit_param = parameter to fit: 0=dg, 1=t0, 2=dg with linear regression
%  ouptut_name = name of output files
%  result_dir = to save results
%  I = use to pre-select inputs (one column, for ALL input files)
%
% Output:
%  r = r-squared between input param and fitted param
% df = number of degrees-of-freedom in fitted model

I = [];
minD = [];
maxN = 2e4;  % max number of input oligos for regression
maxK = 1.2e4; % max number of K-mers for regression
alpha = 0.01;
lasso_exec = 1;
rng('shuffle');


% input data
[idx,D,~,logx,~,bx,tx,k,h] = load_regression_data(data_files,fit_param,I,maxN,alpha);
if (~isempty(minD))
    D((~isnan(tx)).*(tx==min(tx))==1) = minD;
    D((~isnan(D)).*(D<minD)==1) = minD;
end
saveas(h, [result_dir '/input_reporters.jpg'],'jpg');
close all;

% input features
[K,Krows,h] = load_kmer_counts(idx,feature_files,k);
saveas(h, [result_dir '/input_kmers.jpg'],'jpg');
close all;

% FIT: linear regression with L1 regularization
r = 0;
df = 0;
file_name = [result_dir '/run_lasso.out.mat'];
if (exist(file_name, 'file') ~= 2)
    
    j = select_kmers(Krows,K(:,k),D(k),maxK,[result_dir '/kmers_enrichments.mat']);
    if (fit_param == 2)
        [X,Y] = linear_fit_data(K(j,k)',logx(k,:),bx(k),tx(k),tm);
    else
        X = K(j,k)';
        Y = D(k);
    end
    
    % alpha = 1: lasso (L1 reg.), eps: ridge (L2 reg.)
    if (lasso_exec > 0)
        fprintf('running lasso optimization ... (slow)\n');
        [bweights,bfitinfo] = lasso(X,Y,'Alpha',1,'CV',10);
        save(file_name,'bweights','bfitinfo','j','-v7.3');
        fprintf('done\n');
    else
        save([result_dir '/run_lasso.mat'],'X','Y','j','-v7.3');
        system(['rm -rf ' result_dir '/run_lasso.m;' ...
            'echo "load(''run_lasso.mat'');"' ...
            '  > ' result_dir '/run_lasso.m;' ...
            'echo "[bweights,bfitinfo] = lasso(X,Y,''Alpha'',1,''CV'',10);"' ...
            ' >> ' result_dir '/run_lasso.m;' ...
            'echo "save(''run_lasso.out.mat'',''bweights'',''bfitinfo'',''j'',''-v7.3'');"' ...
            ' >> ' result_dir '/run_lasso.m;']);
    end
end

% Feature selection: by regression with L1 regularization
if (exist(file_name, 'file') == 2)
    load(file_name,'bweights','bfitinfo','j');
    x = bfitinfo.Index1SE;
    if (sum(abs(bweights(:,x))>0) < 5)
        x = bfitinfo.IndexMinMSE;
    end
    B0 = bweights(:,x);
    b0 = bfitinfo.Intercept(x);
    i = j;
    i(i==1) = abs(B0) > 0;

    if (fit_param == 2)
        X = linear_fit_data(K(j,k)',logx(k,:),bx(k),tx(k),tm);
    else
        X = K(j,k)';
    end
    
    df = bfitinfo.DF(x);
    fprintf('Lasso: %d features selected\n', df);
    D0 = X*B0 + b0;
    r = corr(D0,D(k));
    fprintf('Lasso: r = %.2f\n', r);
    if (fit_param == 1)
        g = D>1;
        u = unique(D0);
        mint = 0;
        mine = inf;
        for a = 1:max(size(u))
            c = zeros(size(g));
            c(D0>=u(a)) = 1;
            e = sum(c~=g);
            if (e < mine)
                mine = e;
                mint = a;
            end
        end
        g0 = zeros(size(g));
        g0(D0>=u(mint)) = 1;
        fprintf('class: thr = %f, e = %.2f%% (%d,%d)\n',u(mint),100*sum(g0~=g)./size(g,1),sum((g0==1).*(g==0)),sum((g0==0).*(g==1)));
        plot_model_fit_t(D,D0,g,g0,df,[result_dir '/lasso']);
    else
        plot_model_fit_dg(D(k),D0,df,[result_dir '/lasso']);
    end
    
    S = sum(abs(B0));
    F = B0./S;
    write_text_file([result_dir '/lasso_model.txt'], ...
        sortrows([Krows(i) num2cell([B0(abs(B0)>0) F(abs(B0)>0)])],[2,1]));
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    kmer_linear_regression_plot_model(B0(abs(B0)>0),b0,Krows(i),50);
    saveas(h, [result_dir '/lasso_model.jpg'],'jpg');

    save(file_name,'bweights','bfitinfo','j','Krows','i','-v7.3');
    close all;
end



function [idx,D,tm,logx,ax,bx,tx,k,h] = load_regression_data(data_files,fit_param,I,maxN,alpha)

nfiles = max(size(data_files));
ax = [];
bx = [];
tx = [];
rx = [];
px = [];
idx = [];
logx = [];
for j = 1:nfiles
    [a,b,t,~,logY,~,~,~,id,p,r,X] = degradation_model_load(data_files{j});
    tm = X{1}{1};
    p(isnan(b)) = nan;
    
    px = [px; p];
    ax = [ax; a];
    bx = [bx; b];
    tx = [tx; t];
    rx = [rx; r];
    idx = [idx; id];
    if (fit_param == 2)
        logx = [logx; logY];
    else
        logx = [logx; zeros(size(px))];
    end
end

if (fit_param == 1)
    D = tx;
else
    D = -1*ax(:,end);
end

% select a subset for feature selection (lasso)
w = [min(bx,[],2) min(px,[],2) min(rx,[],2)];
i = sum(isnan(bx),2) > 0;
w(i,:) = NaN;
i = max(bx,[],2) < 0.8;
w(i,:) = NaN;
i = max(tx,[],2) > 6;
w(i,:) = NaN;
i = (max(px,[],2) <= 2*alpha).*(max(rx,[],2) <= 0.7) == 1;
w(i,:) = NaN;
if (~isempty(I))
    w(I==0,:) = NaN;
end
w = (w - repmat(prctile(w,1),size(w,1),1));
w = w./repmat(prctile(w,99),size(w,1),1);
w(w<0) = 0;
w(w>1) = 1;

s = sum(w,2);
s = s - min(s);
s = s./max(s);
if (sum(~isnan(s)) > maxN)
    k = s > prctile(s,100*(1-maxN/sum(~isnan(s))));
else
    k = ~isnan(s);
end

% plot regression data
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

subplot(1,2,1);
x = 0:0.01:1;
n = sum(~isnan(w(:,1)));
hold on;
y = hist(w,x);
plot(x,cumsum(y)./n,'linewidth',1.2);
z = hist(s,x);
plot(x,cumsum(z)./n,'linewidth',2);
t = min(s(k));
line([t t],[0 1],'LineStyle','--','color','k','linewidth',1.2);
hold off;
axis tight;
legend({'b' 'p' 'r' 'sum'});
title(sprintf('Reporters: select %d out of %d (%.1f%%)', sum(k),n,100*sum(k)/n));

subplot(1,2,2);
hold on;
if (fit_param == 1)
    x = 1:10;
    y = hist(D(k),x);
    bar(x,y,'facecolor',[0.8 0.8 0.8]);
    x = 1:0.1:10;
    y = hist(D(k),x);
    plot(x,y,'-r','linewidth',1.2);
    xlabel('t0 (hpf)');
else
    x = 0:0.01:1;
    y = hist(D(k),x);
    plot(x,y,'-','linewidth',1.2);
    xlabel('dg rate (1/h)');
end
hold off;
axis tight;
ylabel('number of reporters');


function [K,Krows,h] = load_kmer_counts(idx,feature_files,k)

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
    
    if ((max(size(Kjcol)) ~= max(size(idx))) || (sum(strcmp(Kjcol',idx))<max(size(idx))))
        fprintf('error loading kmer file %s\n', feature_files{j,1});
    else
        K = [K; Kj];
        Krows = [Krows; Kjrows];
    end
end

% plot kmers
h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);

l = cellfun(@length,Krows);
y = full(sum(K(:,k),2))./sum(k);
x = 0:0.1:prctile(y,99.9);
u = unique(l);
hold on;
for i = u'
    z = hist(y(l==i),x);
    plot(x,z./sum(l==i),'-','marker','.','linewidth',1.5);
end
hold off;
axis tight;
xlabel('average occurances (per UTR)');
ylabel('fraction of kmers');
set(gca,'xtick',0:10,'ytick',0:0.2:1,'ylim',[0 1]);
legend(cellstr(strcat('k=', num2str(u))),'box','off');
title(sprintf('Kmers: %d', max(size(Krows))));


function j = select_kmers(Krows,K,D,maxK,enrichment_file_name)

n = size(Krows,1);
if (n > maxK)
    fprintf('selecting top %d of %d kmers by pvalue (%.1f%%)\n', maxK, n, 100*maxK/n);
    if (exist(enrichment_file_name, 'file') ~= 2)
        J = K>0;
        P = kmer_enrichments(J,D, 1,0);
        Q = kmer_enrichments(J,D,-1,0);
        save(enrichment_file_name,'P','Q','Krows');
    else
        load(enrichment_file_name);
    end
    
    v = min([P Q],[],2);
    s = sort(v);
    if (max(size(s)) > maxK)
        s = s(maxK);
    else
        s = s(end);
    end
    j = (v <= s);
else
    j = ones(n,1)>0;
end

fprintf('using %d k-mers (%.1f%%)\n',sum(j),100*sum(j)./n);



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
title(sprintf('df = %d, r = %.2f', df, corr(D0,D)));
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

