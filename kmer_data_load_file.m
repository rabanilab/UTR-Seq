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

function [X, Xrows, Xcols] = kmer_data_load_file(id,filename)

X = [];
Xrows = [];
Xcols = [];
if (exist([filename '.mat'],'file') == 2)
    fprintf('loading kmer file %s\n', filename);
    load([filename '.mat'],'X','Xrows','Xcols');
    if (~isempty(id))
        [~,i] = ismember(id,Xcols);
        Xcols = Xcols(i(i>0));
        X = X(:,i(i>0));
    end
    
elseif (~isempty(id))
    fprintf('reading text file %s\n', filename);
    [X, Xrows, Xcols] = kmer_load_list([filename '.txt'], id);
    save([filename '.mat'],'X','Xrows','Xcols');
end



function [X, Xrows, Xcols] = kmer_load_list(filename,id)

f = fopen(filename);
X = textscan(f,'%s %s %s %s'); %<kmer><count><oligo list><oligo weights>
fclose(f);
X = [X{:}];
Xrows = X(:,1);
Xcols = id';

n = size(X,1);
x = [];
y = [];
v = [];
for k = 1:n
    l = strsplit([X{k,3}],',');
    w = str2double(strsplit([X{k,4}],','));
    [~,j] = ismember(l,Xcols);
    y = [y j];
    i = find(strcmp(Xrows,X{k,1}));
    x = [x i*ones(size(j))];
    v = [v w];
end
X = sparse(x,y,v,size(Xrows,1),size(Xcols,2));
