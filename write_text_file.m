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

function write_text_file(filename, data)
% write text data into a file
% input is a cell

rows = size(data,1);
cdata = to_str(data);

fid = fopen(filename, 'wt');
for (i=1:rows)
    fprintf(fid, '%s\t', cdata{i,1:end-1});
    fprintf(fid, '%s\n', cdata{i,end});
end
fclose(fid);



function h = to_str(v)

[n,m] = size(v);
h = cell(n,m);
for i=1:n
    for j=1:m
        if(isnumeric(v{i,j})+islogical(v{i,j})>0)
            h{i,j} = num2str(v{i,j});
        elseif (isempty(v{i,j}))
            h{i,j} = '';
        else
            h{i,j} = v{i,j};
        end
    end
end
