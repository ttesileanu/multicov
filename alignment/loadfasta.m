function alignment = loadfasta(fname, alphabet, varargin)
% LOADFASTA Load alignment from FASTA file.
%   alignment = LOADFASTA(fname, alphabet) loads an alignment from the
%   FASTA file with the given name, using the provided alphabet. By default,
%   this turns all non-alphabet letters to gaps.
%
%   Options:
%    'keepmask' <v/f>
%       This is either a logical vector showing which positions in the
%       alignment to keep, or a function that is applied to the first
%       alignment sequence to give the mask.
%       (default: keep only upper case letters and '-')
%    'matclean' <b>
%       If this is true, replace invalid characters (for the given alphabet)
%       by gaps. (default: true)
%    'postprocess' <f>
%       A function to apply to the alignment data before calling matclean.
%       If this is empty, the sequences are kept as loaded. (default: [])
%
%   See also ALNMAKE.

% Tiberiu Tesileanu (2012-2014)

% the fastaread version is almost 2x slower and uses more than 50% more
% memory to achieve the same task...

% fasta = fastaread(fname);
% 
% annotations = {fasta.Header};
% 
% alignment = char({fasta.Sequence});
% clear('fasta');
% 
% alignment = matclean(alignment, alphabet);
% alignment = alnmake(alignment, alphabet);
% alignment.annotations = annotations;

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('keepmask', @(c) isstrprop(c, 'upper') | c == '-', ...
    @(p) isvector(p) || isa(p, 'function_handle'));
parser.addParamValue('matclean', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('postprocess', [], @(p) isempty(p) || isa(p, 'function_handle'));

% parse
parser.parse(varargin{:});
params = parser.Results;

% read the file in a cell array of lines
f = fopen(fname);
if f == -1
    error([mfilename ':badfile'], 'Cannot open file.');
end
cdata = textscan(f, '%s', 'delimiter', '\n');
cdata = cdata{1};
fclose(f);

% get the number of sequences
n = sum(cellfun(@(s) ~isempty(s) && s(1) == '>', cdata));

% parse the sequences
data = '';
seq = '';
annotations = cell(1, n);
crt = 0;
for i = 1:length(cdata)
    line = cdata{i};
    if isempty(line)
        continue;
    end
    if line(1) == '>'
        if ~isempty(seq)
            if crt == 1 && isa(params.keepmask, 'function_handle')
                % this is the first sequence: guess the mask
                params.keepmask = params.keepmask(seq);
            end
            seq = seq(params.keepmask);
            if isempty(data)
                data = repmat(blanks(length(seq)), n, 1);
            end
            data(crt, :) = seq;
        end
        
        crt = crt + 1;
        annotations{crt} = line(2:end);
        seq = '';
    else
        % sequences can span several lines
        seq = [seq line]; %#ok<AGROW>
    end
end
% handle the last sequence
if ~isempty(seq)
    if crt == 1 && isa(params.keepmask, 'function_handle')
        % this is the first sequence: guess the mask
        params.keepmask = params.keepmask(seq);
    end
    data(crt, :) = seq(params.keepmask);
end
% make sure we didn't miscount the sequences somehow
data = data(1:crt, :);
annotations = annotations(1:crt);

% do the post processing, if required
if ~isempty(params.postprocess)
    data = params.postprocess(data);
end
% do the cleaning
if params.matclean
    data = matclean(data, alphabet);
end
% make the alignment
alignment = alnmake(data, alphabet);
alignment.annotations = annotations;

end