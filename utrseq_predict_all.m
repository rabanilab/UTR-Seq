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

function utrseq_predict_all(STEP, seq_file_name, bg_seq_file_name, model_file_name, peak_file_name, result_dir)
% This function will use a fitted model, and run its predictions on a test
% set of new sequences (that were not part of the reporter assay)
%
%   STEP 1: predict degradation rates for each input sequence from its
%           kmer counts
%
%   STEP 2: predict short sequence elements that contribution to
%           degradation rate (motifs)
%
% Input Parameters:
%  result_dir = 'example/test_seq_predict_all'
%  seq_file_name = 'example/test_seq.txt'
%  bg_seq_file_name = 'example/bg_seq.txt'
%  model_file_name = 'example/repoters_kmer_regression/run_linear.out.mat'
%  peak_file_name = 'example/repoters_fit_peaks/peaks.mat'
%
% Example:
%   utrseq_predict_all('example/test_seq.txt', 'example/bg_seq.txt', ...
%     'example/repoters_kmer_regression/run_linear.out.mat', ...
%     'example/repoters_fit_peaks/peaks.mat', 'example/test_');

if (exist(result_dir,'dir') ~= 7)
    mkdir(result_dir);
end

% ---------------
% Step 1: predict degradation rates of input sequences (by kmer counts)
% 
% Outputs:
%  1. <result_dir>/deg_rates.mat (predicted degradation rates)
%     <result_dir>/deg_rates.txt
% ---------------
if (STEP > 0)
    utrseq_predict_deg_rates(seq_file_name, bg_seq_file_name, model_file_name, result_dir);
end

% ---------------
% Step 2: predict peak positions in input sequences
% 
% Outputs:
%  1. <result_dir>/peaks.mat (predicted peak positions)
%     <result_dir>/peaks.txt
% ---------------
if (STEP > 1)
    utrseq_predict_peaks(seq_file_name, bg_seq_file_name, model_file_name, peak_file_name, result_dir)
end




function utrseq_predict_deg_rates(seq_file_name, bg_seq_file_name, model_file_name, result_dir)

min_sequence_len = 50;
min_norm_len = 100;
max_norm_len = 120;
minK = 3;
maxK = 7;

%%% Sequences
fprintf('** INPUT sequences\n');
seq_file = [result_dir '/sequences.mat'];
if (~exist(seq_file,'file'))
    [id,seq,seqlen] = load_sequences(seq_file_name,min_sequence_len);
    [~,bgseq] = load_sequences(bg_seq_file_name,min_sequence_len);
    save(seq_file,'id','seq','seqlen','bgseq');
else
    load(seq_file,'id','seqlen');
end
fprintf('Done\n\n');

%%% KMERs
fprintf('** INPUT sequences: KMERs\n');
kmer_dir = regexprep(seq_file_name,'.txt','_kmers');
if (exist(kmer_dir,'dir') ~= 7)
    mkdir(kmer_dir);
end

kmer_files = cellstr(strcat([kmer_dir '/counts.'], num2str((minK:maxK)')));
for i = 1:max(size(kmer_files))
    k = i + minK - 1;
    if (exist([kmer_files{i} '.mat'], 'file') ~= 2)
        extract_kmer_counts(k,seq_file_name,kmer_files{i});
    end
end
fprintf('Done\n\n');

%%% Predict degradation rates
fprintf('** PREDICT: degradation rates\n');
dg_file = [result_dir '/deg_rates.mat'];

if ((min(seqlen) < min_norm_len) + (max(seqlen) > max_norm_len) > 0)
    normalize_length = 1;
    set_bounds = 1;
else
    normalize_length = 0;
    set_bounds = 0;
end

[DG,h] = sequence_to_deg(id,model_file_name,[kmer_dir '/counts'],normalize_length,set_bounds);
saveas(h, [result_dir '/deg_rates.jpg'],'jpg');
save(dg_file,'id','DG');
write_text_file(regexprep(dg_file,'.mat','.txt'),[id num2cell([DG log(2)./DG])]);

close all;
fprintf('Done\n\n');



function utrseq_predict_peaks(seq_file_name, bg_seq_file_name, model_file_name, peak_file_name, result_dir)

min_sequence_len = 50;
minK = 3;
maxK = 7;

%%% Sequences
fprintf('** INPUT sequences\n');
seq_file = [result_dir '/sequences.mat'];
if (~exist(seq_file,'file'))
    [id,seq,seqlen] = load_sequences(seq_file_name,min_sequence_len);
    [~,bgseq] = load_sequences(bg_seq_file_name,min_sequence_len);
    save(seq_file,'id','seq','seqlen','bgseq');
else
    load(seq_file,'id','seq','bgseq');
end
fprintf('Done\n\n');

%%% KMERs
fprintf('** INPUT sequences: KMERs\n');
kmer_dir = regexprep(seq_file_name,'.txt','_kmers');
if (exist(kmer_dir,'dir') ~= 7)
    mkdir(kmer_dir);
end

kmer_files = cellstr(strcat([kmer_dir '/counts.'], num2str((minK:maxK)')));
for i = 1:max(size(kmer_files))
    k = i + minK - 1;
    if (exist([kmer_files{i} '.mat'], 'file') ~= 2)
        extract_kmer_counts(k,seq_file_name,kmer_files{i});
    end
end
fprintf('Done\n\n');

%%% Positional Weights
fprintf('** INPUT sequences: position weights\n');
plot_dir = [result_dir '/weights'];
if (exist(plot_dir,'dir') ~= 7)
    mkdir(plot_dir);
end

wfile = [result_dir '/kmer_weights.mat'];
if (~exist(wfile,'file'))
    W = [];
    S = [];
    for j = 1:max(size(id))
        [Wj,Sj] = kmer_positional_weights(id(j),seq{j},model_file_name);
        W = [W;{Wj}];
        S = [S;Sj];
    end
    save(wfile,'id','W','S');
    
    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    system(['rm -rf ' plot_dir '/*.jpg']);
    ymax = ceil(10*max(abs([W{:}])))/10;
    for j = 1:max(size(id))
        clf;
        plot_weights(seq{j},W{j},regexprep(id{j},'_',' '),S(j),ymax);
        saveas(h, [plot_dir '/' id{j} '.jpg'],'jpg');
    end
    close all;
else
    load(wfile,'W');
end
fprintf('Done\n\n');

%%% Predict peak locations
fprintf('** PREDICT: peak locations\n');
peak_file = [result_dir '/peaks.mat'];

load(peak_file_name,'Mw','Tw','PWM','Q');
PK = [];
A = [];
B = [];
N = [];
for j = 1:max(size(id))
    [PKj,Aid,aj,bj,nj] = assign_peaks(id(j),seq{j},bgseq,W{j},Mw,Tw,Q,PWM);
    PK = [PK;PKj];
    A(j,:) = aj;
    B(j,:) = bj;
    N(j,:) = nj;
end

Xrow = id;
Xid = Aid;
X = sparse(A>0); % genes x motifs
save(peak_file,'Aid','A','B','N','Xid','Xrow','X');
write_text_file(regexprep(peak_file,'.mat','.txt'),PK);

p = cellfun(@isempty,regexp(Xid,':P'))==0;
w1 = full(sum(sum(X(:,p))));
w2 = full(sum(sum(X(:,p==0))));
s = full(sum(X,2));
fprintf('%d positive peaks in %d transcripts (%.0f%%)\n',w1,size(X,1),100*w1/sum(s));
fprintf('%d negative peaks in %d transcripts (%.0f%%)\n',w2,size(X,1),100*w2/sum(s));
fprintf('avg: %.2f peaks per sequence (max=%d)\n',mean(s),max(s));
fprintf('Done\n\n');



function [aid,aseq,alen] = load_sequences(seqfile,min_3utr_len)
% load endogenous sequences

if (nargin < 2)
    min_3utr_len = 50;
end

% load sequences
f = fopen(seqfile);
X = textscan(f,'%s %s');
aid = X{1}; 
aseq = X{2};
fclose(f);
fprintf('LOAD: file = %s\n', seqfile);
fprintf('LOAD: n = %d\n', size(aid,1));

% sort
[~,s] = sort(aid);
aid = aid(s,:);
aseq = aseq(s,:);

% filter by length
alen = cellfun(@length,aseq);
alen(alen<1) = 1;
i = alen>=min_3utr_len;
aid = aid(i,:);
aseq = aseq(i,:);
alen = alen(i,:);
fprintf('LOAD: n = %d after length filter\n', size(aid,1));


function [W,Xid,X,Y,Z] = assign_peaks(id,idseq,bgseq,modelW,Mweight,Tweight,Q,PWM)
% Output:
%  W = list of peaks:
%    [geneid] [position] [peak seq] [peak score] [motif score] [pos/neg] [assigned motif]
%  Xid = motif ids (columns of X)
%  X = abs(peak scores) matrix (ids x motifs)
%    if multiple peaks within one gene, take the sum of their weights
%  Y = peak position matrix (ids x motifs)
%    if multiple peaks within one gene, take the 5'-most
%  Z = peak count matrix (ids x motifs)

Xid = strcat(Q(:,1),':',num2str(cell2mat(Q(:,3))+1));
Xid = regexprep(Xid,':2',':P');
Xid = regexprep(Xid,':0',':N');
Xid = [Xid;'NA:N';'NA:P'];

if (size(modelW,1)>10)
    N = kmer_norm_weights(modelW);
else
    N = modelW;
end
modelW = N./Mweight;

% Z1/Z2: <peak weight> <reporter id> <start> <end> <sequence>
[Z1,~,~,Z2] = peaks_extract(id,idseq,modelW,Tweight);
if (isempty(Z1)+isempty(Z2) == 2)
    W = [];
    X = zeros(size(id,1),size(Xid,1));
    Y = zeros(size(id,1),size(Xid,1));
    Z = zeros(size(id,1),size(Xid,1));
    return;
end

pseq = [];
pscore = [];
pgid = [];
ppos = [];
if (~isempty(Z1))
    pseq = [pseq; Z1(:,5)];
    pscore = [pscore; cell2mat(Z1(:,1))];
    pgid = [pgid; Z1(:,2)];
    ppos = [ppos; cell2mat(Z1(:,3))];
end
if (~isempty(Z2))
    pseq = [pseq; Z2(:,5)];
    pscore = [pscore; cell2mat(Z2(:,1))];
    pgid = [pgid; Z2(:,2)];
    ppos = [ppos; cell2mat(Z2(:,3))];
end

% assign peaks to motifs
[W,Wid] = peaks_assign_seeds(pseq,pscore,bgseq,PWM,Q,0);
W = [pgid num2cell(ppos) W];

% build peak matrix
n = max(size(Xid));
X = zeros(size(id,1),n);
Y = zeros(size(id,1),n);
Z = zeros(size(id,1),n);
if (~isempty(W))
    for i = 1:n
        j = strcmp(Wid,Xid(i));
        if (sum(j)>0)
            [uid,~,t] = unique(W(j,1));
            c1 = accumarray(t,abs(cell2mat(W(j,4))),[],@sum);
            c2 = accumarray(t,cell2mat(W(j,2)),[],@min);
            c3 = accumarray(t,1);
            [j1,j2] = ismember(id,uid);
            X(j1,i) = c1(j2(j2>0));
            Y(j1,i) = c2(j2(j2>0));
            Z(j1,i) = c3(j2(j2>0));
        end
    end
end


function [D1,h,p] = sequence_to_deg(ids, model_file, kmer_file_pref, normbylen, setbounds)
% apply a kmer linear regression model to data
% Input:
%             id = gene ids to run the model on
%     model_file = linear regression model file
% kmer_file_pref = kmers file
%           type = 'dg' or 't0'
%      data_file = fit optimal offset by comparing predicted rates to rates in the data_file

UTRSeq_len = 110;
minL = 1.25*UTRSeq_len;
lowerL = 0.2*UTRSeq_len;

% INPUT: linear regression model
load(model_file,'Krows','B0','b0');
fprintf('Input: model %d df\n', size(Krows,1));

% INPUT: kmer matrix
Klen = cellfun(@length,Krows);
minK = min(Klen);
maxK = max(Klen);
feature_files = cellstr(strcat(kmer_file_pref, '.', num2str((minK:maxK)')));
[X,Xrows] = load_kmer_counts(ids,feature_files);
len = full(sum(X(cellfun(@length,Xrows)==3,:)))' + 2;
len(len<1) = 1;
fprintf('Input: %d kmers, %d ids\n', max(size(Xrows)), max(size(ids)));

% apply model to data
i1 = ismember(Krows,Xrows);
if (sum(i1>0) < size(Krows,1))
    Trows = Krows(i1==0);
    T = zeros(size(Trows,1),size(X,2));
    Xrows = [Xrows;Trows];
    X = [X;T];
end
[~,i2] = ismember(Krows,Xrows); %[Krows(i1) Xrows(i2(i2>0))]
sum(strcmp(Krows,Xrows(i2(i2>0))))
Y = X(i2(i2>0),:)';
D0 = Y*B0 + b0;
D0(len<lowerL) = NaN;

% length correction
if (normbylen)
    [Lz,Dz] = kmer2rate(B0,0,Krows,lowerL,X,Xrows);
    Dz = Dz - floor(prctile(Dz,2));

    h = figure;
    scrsz = get(0,'ScreenSize');
    set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
    xlim = log2([min(Lz) max(Lz)]);
    ylim = [-6 6];
    minD = 0.05;
    k = (~isnan(Dz)).*(Dz>0)==1;
    subplot(2,2,1);
    hold on;
    dscatter(log2(Lz(k)),log2(Dz(k)));
    line(log2([minL minL]),ylim,'LineStyle', '-','color','k');
    line(xlim,log2([minD minD]),'LineStyle', '-','color','k');
    p1 = polyfit(log2(Lz(k)),log2(Dz(k)),1);
    y = polyval(p1,xlim);
    line(xlim,y,'LineStyle', '-','color','r');
    p2 = logistic_fit(log2(Dz(k))',log2(Lz(k))',25,[0 -inf -inf 0],[10 inf inf 20]);
    y = logistic_eval(p2,xlim(1):0.1:xlim(2));
    plot(xlim(1):0.1:xlim(2),y,'-m');
    if (abs(p1(1))<0.2)
        p = [0 0];
    else
        p = p2;
    end
    hold off;
    set(gca,'xlim',xlim,'ylim',ylim);
    axis square;
    xlabel('length');
    ylabel('log dg');
    title(sprintf('raw predictions (n=%d,%.2f)',sum(k),abs(p1(1))));
    subplot(2,2,2);
    hold on;
    if (max(size(p)) == 2)
        z = (Lz./minL).^p(1);
        z(Lz<minL) = 1;
    else
        z = 2.^(logistic_eval(p,log2(Lz)) - logistic_eval(p,0));
    end
    D1 = Dz./z;
    dscatter(log2(Lz(k)),log2(D1(k)));
    line(log2([minL minL]),ylim,'LineStyle', '-','color','k');
    hold off;
    set(gca,'xlim',xlim,'ylim',ylim);
    xlabel('length');
    ylabel('log dg');
    axis square;
    title(sprintf('normalized by length (%d)',max(size(p))));
    subplot(2,2,3);
    x = ylim(1):0.2:ylim(2);
    hold on;
    y = hist(log2(Dz(k)),x);
    plot(y./sum(y),x,'-k','linewidth',2);
    y = hist(log2(D1(k)),x);
    plot(y./sum(y),x,'-r','linewidth',2);
    hold off;
    axis tight;
    ylabel('log dg');
    xlabel('fraction');
    axis square;
    subplot(2,2,4);
    x = xlim(1):0.25:xlim(2);
    hold on;
    y = hist(log2(Lz),x);
    plot(x,y./sum(y),'-k','linewidth',2);
    hold off;
    axis tight;
    xlabel('log2 length');
    ylabel('fraction');
    axis square;

    if (max(size(p)) == 2)
        z = (len./minL).^p(1);
        z(len<minL) = 1;
    else
        z = 2.^(logistic_eval(p,log2(len)) - logistic_eval(p,0));
    end
    D1 = Y*B0;
    D1(len<lowerL) = NaN;
    D1 = D1./z + b0;
else
    h = [];
    p = [];
    D1 = D0;
end

if (setbounds)
    D1 = bounds(D1);
end


function [Lz,Dz] = kmer2rate(B0,b0,Krows,lowerL,Z,Zrows)

Lz = full(sum(Z(cellfun(@length,Zrows)==3,:)))' + 2;
Lz(Lz<1) = 1;
[~,jz] = ismember(Zrows,Krows);
Dz = Z(jz>0,:)'*B0 + b0;
Dz(Lz<lowerL) = NaN;


function [K,Krows] = load_kmer_counts(idx,feature_files)

K = [];
Krows = [];
for j = 1:max(size(feature_files))
    [Xj, Xjrows, Kjcol] = kmer_data_load_file([],feature_files{j});
    if (isempty(Xj))
        fprintf('error loading kmer file %s\n', feature_files{j});
        continue;
    end
    
    if (~isempty(idx))
        [~,i] = ismember(idx,Kjcol);
        Xi = zeros(size(Xj,1),size(i,1));
        Xi(:,i>0) = Xj(:,i(i>0));
        K = [K; Xi];
    else
        K = [K; Xj];
    end
    Xjrows = regexprep(Xjrows,'_','.');
    Krows = [Krows; Xjrows];
end

function d0 = bounds(D0)

d0 = -1*degradation_model_set_bounds(-1*D0);


function [X,Xrows] = kmer_counts(seq,k)

if (iscell(seq))
    m = max(size(seq));
    Xrows = cellstr(kmer_all(k));
    X = zeros(max(size(Xrows)),m);
    for j = 1:m
        X(:,j) = kmer_counts_str(seq{j},k);
    end
else
    [X,Xrows] = kmer_counts_str(seq,k);
end

function [X,Xrows] = kmer_counts_str(seq,k)
[m,n] = size(seq);

% split sequence into kmers of size k
K = strings(n-k+1,m);
for j = 1:m
    T = repmat(char(0),k,n-k+1);
    for i = 1:k
        T(i,:) = seq(j,1+(i-1):n-(k-i));
    end
    K(1:n-k+1,j) = cellstr(T');
end

% kmer counts
Xrows = cellstr(kmer_all(k));
X = zeros(max(size(Xrows)),m);
for j = 1:m
    [u,~,t] = unique(K(:,j));
    n = accumarray(t,1);
    [i1,i2] = ismember(Xrows,u);
    X(i1,j) = n(i2(i2>0));
end


function plot_weights(seq,w,id,dg,ymax)

c = [0.00 0.45 0.74; ... %b
     0.20 0.75 0.93; ... %b
     0.47 0.67 0.19; ... %g
     0.00 0.45 0.05; ... %g
     0.64 0.08 0.18; ... %r
     0.85 0.13 0.10; ... %r
     0.93 0.69 0.25; ... %y
     0.90 0.85 0.03];    %y
unique(colormap('lines'),'rows');
c = [c;c];

% t = 0.05;
% if (sum(w>t)>0)
%     [cellstr(repmat(id,sum(w>t),1)) num2cell([find(w>t)' w(w>t)']) cellstr(seq(w>t)')]
% end

[m,n] = size(w);
hold on;
L = cell(1,m);
for i = 1:m
    plot(1:n,w(i,:),'-','color',c(i,:),'linewidth',1.2);
    L{i} = sprintf('(dg = %.2f, %.1f min)', dg(i),60*log(2)/dg(i));
end
line([1 n],[0 0],'LineStyle', '-','color','k');
for i = 1:n
    fs = round(1400/n);
    fs(fs<2) = 2;
    text(i,0,seq(i),'FontSize',fs,'FontName','Courier');
end
hold off;
axis tight;
set(gca,'ylim',[-1*ymax ymax]);
xlabel('position');
ylabel('kmer weight');
title(sprintf('%s (%d nt)',id,n));
legend(L,'box','off');
box off;

