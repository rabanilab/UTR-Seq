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
% This is an example of how to run the entire pipeline to analyze
% UTR-Seq results
% ----------------------------------------------------------------------

function utrseq_fit_all(STEP, seq_file_name, data_file_name, tm, result_base)
% ----------------------------------------------------------------------
% This function will fit a linear regression model to predict short
% regulatory sequences from a massively parallel reporter assay (MPRA) of
% mRNA stability. There are 3 steps:
%
%   STEP 1: fit an exponential decay degradation rate for each reporter
%           from its RNA-Seq counts
%
%   STEP 2: fit a linear regression model that converts k-mer counts in reporter
%           sequences into degradation rates
%
%   STEP 3: identify short sequences with maximal contribution to
%           degradation rate (peaks), and combine similar sequences into
%           PWM models (motifs).
%
% Input Parameters:
%  result_base = 'example/repoter_'
%  seq_file_name = 'example/reporter_seq.txt'
%  data_file_name = 'example/reporter_counts.txt'
%  tm = [0 1 2 3 4 5 6 7 8 10 12];
%
% Example:
%  utrseq_fit_all(3, 'example/reporter_seq.txt', 'example/reporter_counts.txt',
%    [0 1 2 3 4 5 6 7 8 10 12], 'example/reporter_');
% ----------------------------------------------------------------------


% ---------------
% Step 1: fit degradation rates to temporal samples
%
% Outputs:
%  1. <result_base>_fit_degradation/degradation.norm.mat (normalized counts)
%     <result_base>_fit_degradation/degradation.norm.txt
%  2. <result_base>_fit_degradation/degradation.fit.mat (fitted degradation model)
%     <result_base>_fit_degradation/degradation.fit.txt
% ---------------
if (STEP > 0)
    output_dir = [result_base 'fit_degradation'];
    run_fit_degradation_model(tm, data_file_name, output_dir);
end

% ---------------
% Step 2: fit kmer regression model to sequences and degradation rates
%
% Outputs:
%  1. <result_dir>/reporter_seq_kmers (kmer counts per reporters)
%  1. <result_base>_kmer_regression/run_lasso.out.mat (lasso fit regression)
%  2. <result_base>_kmer_regression/run_linear.out.mat (linear fit regression model)
% ---------------
if (STEP > 1)
    input_file = [result_base '_fit_degradation/degradation.fit.mat'];
    output_dir = [result_base '_kmer_regression'];
    run_fit_kmer_regression(input_file, seq_file_name, output_dir);
end

% ---------------
% Step 3: compute peaks
%
% Outputs:
%  1. <result_base>_fit_peaks/peaks.mat (peak PWM models)
%  2. <result_base>_fit_peaks/peaks.logo.* (peak PWM logos)
%  3. <result_base>_fit_peaks/peaks.a.txt (all peaks, with scores)
% ---------------
if (STEP > 2)
    model_file = [result_base '_kmer_regression/run_linear.out.mat'];
    output_dir = [result_base '_fit_peaks'];
    run_fit_regression_peaks(model_file, seq_file_name, output_dir);
end
