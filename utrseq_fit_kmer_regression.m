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

function utrseq_fit_kmer_regression(input_file, input_seq_file, result_dir)

%%% kmer data
kmer_dir = regexprep(input_seq_file,'.txt','_kmers');
if (exist(kmer_dir,'dir') ~= 7)
    mkdir(kmer_dir);
end

minK = 3;
maxK = 7;
kmer_files = cellstr(strcat([kmer_dir '/counts.'], num2str((minK:maxK)')));

load(input_file,'id');
for i = 1:max(size(kmer_files))
    k = i + minK - 1;
    if (exist([kmer_files{i} '.txt'], 'file') ~= 2)
        system(['./extract_kmer_counts.pl ' input_seq_file ' -k ' num2str(k) ' > ' kmer_files{i} '.txt;']);
    end        
    if (exist([kmer_files{i} '.mat'], 'file') ~= 2)
        kmer_data_load_file(id,kmer_files{i});
    end
end

%%% regression
if (exist(result_dir,'dir') ~= 7)
    mkdir(result_dir);
end

% feature selection & LASSO optimization
kmer_lasso_regression_cont({input_file},kmer_files,0,result_dir);

% linear fit
lasso_file = [result_dir '/run_lasso.out.mat'];
kmer_linear_regression_cont({input_file},kmer_files,lasso_file,0,result_dir);
