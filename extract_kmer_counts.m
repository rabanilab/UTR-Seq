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

function extract_kmer_counts(K,sequence_file,filename)

f = fopen(sequence_file);
X = textscan(f,'%s %s');
fclose(f);
id = X{1};
seq = char(upper(X{2}));
seq = cellstr(seq);

% count k-mers
Xrows = all_kmers(K);
Xcols = id;
X = sparse(size(Xrows,1),size(Xcols,1));
for i = 1:size(seq,1)
    ci = nmercount(seq{i},K,0);
    if (isempty(Xrows))
        Xrows = ci(:,1);
        X = [X cell2mat(ci(:,2))];
    else
        [x,y] = ismember(Xrows,ci(:,1)); % [Xrows(x) ci(y(y>0),1)]
        X(x,i) = cell2mat(ci(y(y>0),2));
    end
end

% remove empty k-mers
i = sum(X,2) > 0;
fprintf('size %d K-mers (%d sequences), %d K-mers exist\n', K, size(id,1), full(sum(i)));
X = X(i,:);
Xrows = Xrows(i);

save([filename '.mat'],'X','Xrows','Xcols');

function S = all_kmers(K)

L = ['A'; 'C'; 'G'; 'T'];

S = [];
if (K == 1)
    S = cellstr(L);
elseif (K > 1)
    Si = all_kmers(K-1);
    for i = 1:max(size(L))
        S = [S;strcat(L(i),Si)];
    end
end
    