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

function utrseq_fit_peaks(model_file, sequence_file, result_dir)
% compute regression peaks using sequences in the input sequence file

if (exist(result_dir,'dir') ~= 7)
    mkdir(result_dir);
end

% load sequences
fprintf('loading sequences ... \n');
seq_file = [result_dir '/sequences.mat'];
if (exist(seq_file, 'file') ~= 2)
    [ids, sequences] = load_sequences(sequence_file);
    sequences = char(sequences);
    save(seq_file,'ids','sequences');
else
    load(seq_file,'ids','sequences');
end
fprintf('Done\n\n');

% positional weights
weights_file = [result_dir '/weights.mat'];
if (exist(weights_file,'file') ~= 2)
    fprintf('calculating sequence weights ... \n');
    W = kmer_positional_weights(ids,sequences,model_file);
    save(weights_file,'ids','W');
    fprintf('Done\n\n');
else
    load(weights_file,'W');
end

% find regression peaks
fprintf('searching regression peaks ... \n');
peak_file = [result_dir '/peaks.mat'];
if (exist(peak_file, 'file') ~= 2)
    [Z,Zu,Zw,Q,PWM,Mw,Tw] = fit_peaks(ids,sequences,W,result_dir);
    save(peak_file,'Z','Zu','Zw','Q','PWM','Mw','Tw');
else
    load(peak_file,'Q','PWM');
end
fprintf('Done\n\n');

% peak logos
delete([result_dir '/peaks.logo.*']);
for j = 1:max(size(PWM))
    fprintf('%s\n',Q{j,1});
    if (Q{j,3}>0)
        fname = ['P.' Q{j,1}];
    elseif (Q{j,3}<0)
        fname = ['N.' Q{j,1}];
    end
    write_text_file([result_dir '/peaks.logo.' fname '.txt'], ...
        [{'A';'C';'G';'T'} num2cell(PWM{j})]);
    
    h = pwm_logo(PWM{j});
    saveas(h,[result_dir '/peaks.logo.' fname '.jpg'],'jpg');
    close(h);
end


function [id,seq] = load_sequences(sequence_file,IDs)

if (nargin < 2)
    IDs = [];
end

f = fopen(sequence_file);
X = textscan(f,'%s %s');
fclose(f);
id = X{1};
seq = char(upper(X{2}));
seq = cellstr(seq);

if (~isempty(IDs))
    [i,j] = ismember(IDs,id);
    S = cell(size(IDs));
    S(i) = seq(j(j>0));
    if (sum(i==0) > 0)
        S{i==0} = '';
    end
    id = IDs;
    seq = S;
end


function [Z,Zu,Zw,Q,PWM,Mweight,Tweight] = fit_peaks(id,idseq,W,result_dir)
% Q = [consensus] [number of peaks] [POS/NEG]
% PWM = position weight matrix per peak motif

iterative_fit = 1;
peak_prctile = 99;
maxf = 0.01;

% positional weights
Wn = kmer_norm_weights(W);
Mweight = max(max(abs(Wn)));
wnn = Wn./Mweight;

% calculate thresholds for peak selection
xj = abs(wnn(:,sum(abs(wnn),1)>0));
dw = 5;
d1 = xj(:);
t1 = prctile(d1,peak_prctile);
d2 = [];
for j = 1:(size(xj,2)-dw+1)
    d2 = [d2 mean(xj(:,j:(j+dw-1)),2)];
end
d2 = d2(:);
t2 = prctile(d2,peak_prctile);
Tweight = t2;
fprintf('fit peaks: pct = %.1f, thr = %.2f, mx = %.2f\n', 100-peak_prctile, Tweight, Mweight);

% extract peak sequences
[~,~,~,~,~,~,Z,Zu,Zw] = peaks_extract(id,idseq,wnn,Tweight);
write_text_file([result_dir '/peaks.a.txt'],Z);
write_text_file([result_dir '/peaks.u.txt'],Zu);
write_text_file([result_dir '/peaks.w.txt'],Zw);

% motifs in peak sequences
pseq = Zw(:,2);
pscore = cell2mat(Zw(:,1));
Q = [];
PWM = [];

% positive peaks
k = pscore > 0; 
fprintf('model: %d positive peaks\n', sum(k));
if (sum(k)>0)
    [Qpos,PWMpos,~,L1] = select_seeds(pseq(k),pscore(k),4,6);
    Q = [Q;[Qpos(:,[2 5]) num2cell(ones(size(Qpos,1),1))]];
else
    L1 = [];
    Qpos = [];
    PWMpos = [];
end

% negative peaks
k = pscore < 0; 
fprintf('model: %d negative peaks\n', sum(k));
if (sum(k)>0)
    [Qneg,PWMneg,~,L2] = select_seeds(pseq(k),pscore(k),4,6);
    Q = [Q;[Qneg(:,[2 5]) num2cell(-1*ones(size(Qneg,1),1))]];
else
    L2 = [];
    Qneg = [];
    PWMneg = [];
end

% select peaks by weight
PWM = [PWMpos;PWMneg];
L = [L1;L2];
w = cell2mat(L(:,1));
[u,~,t] = unique(strcat(L(:,3),':',cellstr(num2str(w>0))));
d = cellfun(@isempty,regexp(u,':0','once'));
s = accumarray(t,abs(w));

minW = maxf*sum(s);
k = zeros(size(Q,1),1);
i = find(s>=minW);
npos = 0;
nneg = 0;
for j = 1:max(size(i))
    uid = strtok(u(i(j)),':');
    k = k + strcmp(Q(:,1),uid);
    npos = npos + (d(i(j))==1);
    nneg = nneg + (d(i(j))==0);
end
if (npos == 0)
    [~,i] = max(s.*(d==1));
    uid = strtok(u(i),':');
    k = k + strcmp(Q(:,1),uid);
    npos = 1;
end
if (nneg == 0)
    [~,i] = max(s.*(d==0));
    uid = strtok(u(i),':');
    k = k + strcmp(Q(:,1),uid);
    nneg = 1;
end
fprintf('SELECT: %d positive, %d negative (w >= %.2f)\n',npos,nneg,minW);
Q(k>0,:)
    
% reassign peaks
[W0,Wid0,xPWM0,xPWMstat0,xPWMw0] = peaks_assign_seeds(pseq,pscore,pseq,PWM(k>0),Q(k>0,:));
[xPWMstat0 num2cell(xPWMw0)]

if (iterative_fit)
    i1 = cell2mat(xPWMstat0(:,3)) == 1;
    j1 = cell2mat(W0(:,4)) == 1;
    [PWM1,Q1,S1] = fit_seeds(pseq(j1,:),pscore(j1,:),W0(j1,:),Wid0(j1),xPWM0(i1),xPWMstat0(i1,:),xPWMw0(i1,:),minW);
    [S1,i] = sort(S1);
    PWM1 = PWM1(i);
    Q1 = Q1(i,:);
    [Q1 num2cell(S1)]
    
    c = 1;
    n = size(Q1,1);
    while (c < n)
        i1 = ones(size(Q1,1),1);
        i1(c) = 0;
        i1 = (i1==1);
        [W1,Wid1,xPWM1,xPWMstat1,xPWMw1] = peaks_assign_seeds(pseq(j1,:),pscore(j1,:),pseq(j1,:),PWM1(i1),Q1(i1,:));
        [PWM1,Q1,S1] = fit_seeds(pseq(j1,:),pscore(j1,:),W1,Wid1,xPWM1,xPWMstat1,xPWMw1,minW);
        if (size(Q1,1)>=n)
            c = c+1;
        else
            n = size(Q1,1);
        end
    end
    [Q1 num2cell(S1)]
    
    i2 = cell2mat(xPWMstat0(:,3)) == -1;
    j2 = cell2mat(W0(:,4)) == -1;
    [PWM2,Q2,S2] = fit_seeds(pseq(j2,:),pscore(j2,:),W0(j2,:),Wid0(j2),xPWM0(i2),xPWMstat0(i2,:),xPWMw0(i2,:),minW);
    [S2,i] = sort(S2);
    PWM2 = PWM2(i);
    Q2 = Q2(i,:);
    [Q2 num2cell(S2)]

    c = 1;
    n = size(Q2,1);
    while (c < n)
        i2 = ones(size(Q2,1),1);
        i2(c) = 0;
        i2 = (i2==1);
        [W2,Wid2,xPWM2,xPWMstat2,xPWMw2] = peaks_assign_seeds(pseq(j2,:),pscore(j2,:),pseq(j2,:),PWM2(i2),Q2(i2,:));
        [PWM2,Q2,S2] = fit_seeds(pseq(j2,:),pscore(j2,:),W2,Wid2,xPWM2,xPWMstat2,xPWMw2,minW);
        if (size(Q2,1)>=n)
            c = c+1;
        else
            n = size(Q2,1);
        end
    end
    [Q2 num2cell(S2)]
    
    PWM = [PWM1;PWM2];
    Q = [Q1;Q2];
else
    PWM = xPWM0;
    Q = xPWMstat0;
end

if (~isempty(result_dir))
    write_text_file([result_dir '/peaks.kmers.txt'],[L1;L2]);
end

function [pwmC,pseudoC,weightN,minI,maxIter,probWeight] = opt_param()

pwmC = 20;
pseudoC = 1;
weightN = 1e3;
minI = 0.1;
maxIter = 25;
probWeight = 0.1;

function [PWM,Q,W] = fit_seeds(pseq,pscore,W0,Wid0,xPWM0,xPWMstat0,xPWMw0,minW)

[pwmC,pseudoC,weightN,minI,maxIter] = opt_param();

S0 = -inf;
S1 = sum(cell2mat(W0(:,3)));
iter = 0;
while ((S1>S0)*(iter<maxIter) == 1)
    PWM = xPWM0;
    Q = xPWMstat0;
    W = xPWMw0;
    fprintf('[iter %d] %.2e\n', iter, S1);

    Xid = unique(Wid0);
    xPWMid = strcat(xPWMstat0(:,1),':',num2str(cell2mat(xPWMstat0(:,end))+1));
    xPWMid = regexprep(xPWMid,':2',':P');
    xPWMid = regexprep(xPWMid,':0',':N');
    n = size(Xid,1);
    xPWM1 = cell(n,1);
    xPWMstat1 = cell(n,3);
    xPWMw = zeros(n,1);
    for i = 1:n
        j1 = strcmp(Wid0,Xid{i})==1;
        j2 = strcmp(xPWMid,Xid{i})==1;
        if (sum(j2)>0)
            [~,~,sseq] = peaks_pwm_score(W0(j1,1),xPWM0{j2},pwmC,pseudoC);
            [xPWM1{i},~,cons] = pwm_construct(W0(j1,1),abs(cell2mat(W0(j1,2))),[],cell2mat(sseq),weightN,minI);
        else
            [xPWM1{i},~,cons] = pwm_construct(W0(j1,1),abs(cell2mat(W0(j1,2))),[],[],weightN,minI);
        end
        xPWMstat1(i,:) = [cons num2cell([sum(j1) 2*isempty(regexp(Xid{i},':N','once'))-1])];
        xPWMw(i) = sum(sum(abs(cell2mat(W0(j1,2)))));
        fprintf('[%d] %s -> %s (w = %.2f)\n', i,Xid{i},cons,xPWMw(i));
    end
    
    q = (xPWMw>=minW);
    if (sum(q) == 0)
        [~,q] = max(xPWMw);
    end

    [W0,Wid0,xPWM0,xPWMstat0,xPWMw0] = peaks_assign_seeds(pseq,pscore,pseq,xPWM1(q),xPWMstat1(q,:));
    S0 = S1;
    S1 = sum(cell2mat(W0(:,3)));
    iter = iter + 1;
end

function [Q,PWM,A,seq_list] = select_seeds(pseq,pscore,mink,maxk)
% select seeds by score
% Q = [seed] [consensus] [seed length] [seed score] [number of peaks]

[~,~,weightN,minI] = opt_param();
min_opt = 100;

% count kmers in peak sequences
Cid = [];
C = [];
plen = cellfun(@length,pseq);
for j = mink:maxk
    [kid,kc] = kmer_counts(pseq,j);
    Cid = [Cid;kid];
    C = [C;kc];
end

% weight[kmer,peak] = (fraction of kmer length out of overall peak length)*(peak score)
W = cellfun(@length,Cid)*abs(pscore./plen)';

% select seeds and optimize motifs
Kseed = [];
Kscore = [];
Kcnt = [];
Kcons = [];
PWM = [];
A = [];
seq_list = [];

% kmer score = sum of weights of all peaks that contain a kmer,
% normalized by kmer length relative to peak's length
I = (ones(1,size(C,2)) == 1);
kmerScore = sum((C(:,I)>0).*W(:,I),2);

while ((sum(I)>0)*(sum(kmerScore>0)>0) > 0)
        
    % select a kmer with maximal score
    [s,j] = max(kmerScore);
    q = (C(j,:)>0).*(I>0) == 1;
    W(j,:) = 0;
    i = (sum((C>0)==repmat(q,size(C,1),1),2) == size(C,2));
    W(i,:) = 0;
    fprintf('seed = %s (score = %.2e) %d\n', Cid{j},s,sum(q));
    
    % calculate a consensus sequence
    [pwm,align,cons] = pwm_construct(pseq(q),abs(pscore(q)),Cid(j),[],weightN,minI);    
%     if (sum(I) >= min_opt)
%         [pwm,align,cons,q] = extend_cons(pwm,cons,q,pseq,pscore,I);
%     end
    
    Kseed = [Kseed;Cid(j)];
    Kscore = [Kscore;s];
    Kcons = [Kcons;{cons}];
    Kcnt = [Kcnt;sum(q>0)];
    PWM = [PWM;{pwm}];
    A = [A;{align}];
    seq_list = [seq_list;[num2cell(pscore(q)) pseq(q) repmat(Kcons(end),sum(q),1) repmat(Cid(j),sum(q),1)]];
        
    % keep only peaks that were not selected
    I = I.*(q==0) == 1;
    kmerScore = sum((C(:,I)>0).*W(:,I),2);
end
if (sum(I)>0)
    seq_list = [seq_list;[num2cell(pscore(I)) pseq(I) repmat({'NA'},sum(I),1) repmat({'NA'},sum(I),1)]];
end

Q = [Kseed Kcons num2cell([cellfun(@length,Kseed) Kscore Kcnt])];
[Q,i] = sortrows(Q,[-5 -4 -3]);
if (~isempty(PWM))
    PWM = PWM(i);
    A = A(i);
end
fprintf('selected seeds: %d\n', size(Q,1))
j = min(size(Q,1),20);
n = [Q{:,end}]';
fprintf('top %d seeds:\n', j)
[Q(1:j,:) num2cell(round(100*n(1:j)./sum(n)))]

function [pwm,align,cons,q] = extend_cons(pwm,cons,q,pseq,pscore,I)

if (nargin < 6)
    I = ones(size(pscore))' == 1;
end

[pwmC,pseudoC,weightN,minI,maxIter,prob_weight] = opt_param();

iter = 0;
i1 = q;
i0 = zeros(size(i1));
while ((sum(i1~=i0) > 0)*(iter <= maxIter) == 1)
    [prob,~,sseq] = peaks_pwm_score(pseq,pwm,pwmC,pseudoC);
    prob = log2(prob);
    min_prob = select_opt_threshold(prob,(I+i1)',prob_weight);
    fprintf('[%d] %s (p > %.2e) %d\n', iter,cons,min_prob,sum(i1));
    i0 = i1;
    [pwm,align,cons] = pwm_construct(pseq(i0),abs(pscore(i0)),[],cell2mat(sseq(i0)),weightN,minI);
    i1 = (prob' > min_prob).*(I>0) + q > 0;
    iter = iter + 1;
end
fprintf('[L] %s (p > %.2e) %d\n\n', cons,min_prob,sum(i1));
q = i0;

function [T,mine] = select_opt_threshold(S,I,prob_weight)

[S,i] = sort(S);
q2 = cumsum(I(i)==2); % sum 1s <=T
q1 = sum(I(i)==1) - cumsum(I(i)==1); % sum 1s > T
q0 = sum(I(i)==0) - cumsum(I(i)==0); % sum 0s > T
e = q2./(sum(q2)+1) + q1./(sum(q1)+1) + prob_weight*q0./(sum(q0)+1);
[mine,j] = min(e);
T = S(j);

% clf;
% hold on;
% plot(q2./sum(q2),'-r','linewidth',1.5);
% plot(q1./sum(q1),'-c','linewidth',1.5);
% plot(q0./sum(q0),'-b','linewidth',1.5);
% plot(e,'-k','linewidth',1.5);
% hold off;
% axis tight;

% clf;
% x = -30:2:max(S)+5;
% y2 = hist(S(I(i)==2),x);
% y2 = y2./sum(y2);
% y1 = hist(S(I(i)==1),x);
% y1 = y1./sum(y1);
% y0 = hist(S(I(i)==0),x);
% y0 = y0./sum(y0);
% hold on;
% plot(x,y2,'-r','linewidth',1.5);
% plot(x,y1,'-c','linewidth',1.5);
% plot(x,y0,'-b','linewidth',1.5);
% line([T T],[0 max([y2 y1 y0])],'LineStyle','-','color','k','linewidth',1.5);
% hold off;
% xlabel('log P');
% legend({'pos' 'unknown' 'others'},'location','bestOutside');
% axis tight;

function [Cid,C] = kmer_counts(S,k)
% Count overlapping k-mers in input sequences

n = size(S,1);
Kseq = [];
Kid = [];
for j = 1:n
    s = S{j};
    l = length(s);
    for i = 1:l-k+1
        Kseq = [Kseq; s(i:i+k-1)];
        Kid = [Kid;j];
    end
end

if (~isempty(Kseq))
    Kseq = cellstr(Kseq);
    [Cid,~,t1] = unique(Kseq);
    W(:,1) = t1;
    W(:,2) = Kid;
    [W,~,t2] = unique(W,'rows');
    c = accumarray(t2,1);
    C = zeros(max(W(:,1)),n);
    for i = 1:size(c,1)
        C(W(i,1),W(i,2)) = c(i);
    end
else
    C = [];
    Cid = [];
end

function [N,A,C] = pwm_construct(kmers,weights,seed,A,weightN,minI)
% generate a PWM matrix from input kmers by weights
%
% Input:
%  kmers = kmer sequences
%  weights = kmer weights
%  seed = alignment seed (if exists)
%  weightN = weighted profile normalization factor
%  minI = minimal position information content to keep in final PWM
% Output:
%  N = PWM matrix
%  A = multiple sequence alignment
%  C = consensus sequence

if (nargin < 2)
    weights = ones(size(kmers));
end
if (nargin < 3)
    seed = [];
end
if (nargin < 4)
    A = [];
end
if (nargin < 5)
    weightN = 1e3;
end
if (nargin < 6)
    minI = 0.1;
end


% multiple sequence alignment
if (isempty(A))
    if (size(kmers,1)==2)
        A = multialign([kmers;kmers],'terminalGapAdjust',true,'GapOpen',100);
        A = A(1:2,:);
    elseif (~isempty(seed))
        A = seedalign(kmers,seed);
    elseif (size(kmers,1)>2)
        A = multialign(kmers,'terminalGapAdjust',true,'GapOpen',100);
    else
        A = kmers{1};
    end
end
[n,m] = size(A);

% weighted counts per position, PWM
L = {'A' 'C' 'G' 'T' '-'};
N = zeros(5,m);
for j = 1:m
    [u,~,t] = unique(A(:,j));
    c = accumarray(t,weights);
    c = round(c*weightN);
    c(c<1) = 1;
    for k = 1:max(size(u))
        f = (strcmp(L,u(k)));
        N(f,j) = c(k);
    end
end
N = N(1:4,:) + repmat(round(N(end,:)/4),4,1);
N = N./repmat(sum(N,1),4,1);

% calculate information content
I = seqlogo(N,'DISPLAYLOGO','FALSE');
I = sum(I{2},1)';
k = find(I>minI,1,'first'):find(I>minI,1,'last');
N = N(:,k);
A = A(:,k);

% consensus
C = '';
for j = 1:size(N,2)
    [~,k] = max(N(:,j));
    C = [C L(k)];
end
C = [C{:}];

function A = seedalign(kmers,seed)

l = cellfun(@length,kmers);
f = cellfun(@min,strfind(kmers,seed));
s = max(f)-f;
e = l+s;
e = max(e)-e;

n = size(kmers,1);
for i = 1:n
    A(i,:) = [repmat('-',1,s(i)) kmers{i,:} repmat('-',1,e(i))];
end

function h = pwm_logo(PWM,N,AB)

if (nargin < 2)
    N = 1000;
end
if (nargin < 3)
    AB = 'ACGU';
end

% generate random sequences
rng('shuffle');
S = [];
for i = 1:size(PWM,2)
    c = cumsum(PWM(:,i));
    r = rand(N,1);
    p = zeros(size(r));
    for j = size(PWM,1):-1:1
        p(r<=c(j)) = j;
    end
    S = [S AB(p)'];
end

% display logo
[~,h] = seqlogo(S);
