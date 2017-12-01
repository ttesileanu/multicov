function [newalign, alignfrac, aligndata] = seqtruncate(alignment, varargin)
% SEQTRUNCATE Truncate and/or find mapping between alignment columns and
% reference sequence.
%   newalign = SEQTRUNCATE(alignment, idx) truncates the alignment to the
%   positions where the idxth sequence is not gapped, and updates the
%   refseq member of the alignment to correspond to these positions.
%
%   newalign = SEQTRUNCATE(alignment) uses the first sequence as the
%   reference.
%
%   SEQTRUNCATE(..., 'seq', seq) uses the given sequence as reference
%   sequence, and performs an alignment to the idxth sequence to find the
%   mapping between alignment columns and this reference sequence. Only the
%   positions for which a match exists are kept. 'Seq' can be a cell array
%   of sequences, one for each alphabet in the alignment.
%
%   SEQTRUNCATE(..., 'truncate', false) performs the match but does not do
%   the truncation, leaving those positions that do not have a non-gap
%   match to the reference sequence marked as 0.
%
%   [newalign, alignfrac] = SEQTRUNCATE(...) also returns alignfrac, the
%   fraction of residues that appear in both the alignment and the
%   reference sequence that agree between the two. This is 1 if an exact
%   match was found.
%
%   [newalign, alignfrac, aligndata] = SEQTRUNCATE(...) returns as the
%   third output a cell array of alignment results, in nwalign format, one
%   element for each alphabet.
%
%   Options:
%    'posnames' <v/c>
%       Array of identifiers for each position in the reference sequence.
%       This can be either a single array (for single-alphabet alignments),
%       or a cell array of arrays, one cell per alphabet (which works for
%       alignments with 1 or more alphabets). The identifiers can be
%       numeric or strings. If one of the cells is empty, the default
%       identifiers are used.
%       (default: numeric indices starting at 1)
%    'seq' <s/c>
%       A string or a cell array of strings (when the alignment has
%       multiple alphabets) giving the reference sequence to be used.
%       (default: use the idxth or the first sequence in the alignment)
%    'truncate' <b>
%       Whether to perform the truncation of the alignment to columns where
%       there is a (non-gap) match to the reference sequence, or to just
%       mark these as 0 in refseq.
%       (default: true)
%    'verbose' <b>
%       If this is true, the function warns if the sequence found in the
%       alignment is not a perfect match to the reference sequence. This
%       only works if 'seq' is provided.
%       (default: true)
%
% See also: SEQSEARCH.

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('idx', 1, @(i) isscalar(i) && isnumeric(i) && isreal(i) && i == floor(i));

parser.addParamValue('seq', [], @(s) (ischar(s) && isvector(s)) || (iscell(s) && all(cellfun(@(ss) ischar(ss) && isvector(ss), s))));
parser.addParamValue('truncate', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('posnames', {}, @(c) isvector(c) && (iscell(c) || isnumeric(c)));
parser.addParamValue('verbose', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% find out where each alphabet starts and where it ends -- we'll need this later
ranges = getalpharanges(alignment);

% for later use: find out what the gap characters are for each of the
% alignments
gapchars = blanks(length(alignment.alphabets));
for i = 1:length(alignment.alphabets)
    alphabet = alignment.alphabets{i};
    if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
        alphabet = alphabet(4:end);
    end
    letters = alphagetletters(alphabet);
    gapchars(i) = letters(1);
end    

% check whether we were given a sequence
if ~isempty(params.seq) && ~iscell(params.seq)
    params.seq = {params.seq};
elseif isempty(params.seq)
    % if not, use one from the alignment
    params.seq = cell(1, size(ranges, 2));
    for i = 1:size(ranges, 2)
        subseq = alignment.data(params.idx, ranges(1, i):ranges(2, i));
        gapch = gapchars(i);
        params.seq{i} = subseq(subseq ~= gapch);
    end
end

newalign = alignment;

% transform alphabets to a form that nwalign likes
nalphas = length(alignment.alphabets);
swalphas = cell(1, nalphas);
for i = 1:nalphas
    alphabet = alignment.alphabets{i};
    if length(alphabet) >= 3 && strcmp(alphabet(1:3), 'gap')
        alphabet = alphabet(4:end);
    end
    if strcmp(alphabet, 'protein')
        swalphas{i} = 'AA';
    elseif strcmp(alphabet, 'dna') || strcmp(alphabet, 'rna')
        swalphas{i} = 'NT';
    else
        error([mfilename ':badalpha'], ['Unsupported alphabet ''' alphabet '''.']);
    end
end

% find the match to the subsequence for each alphabet
aligndata = cell(1, nalphas);
for j = 1:nalphas
    subseq = alignment.data(params.idx, ranges(1, j):ranges(2, j));
    % need the ungapped sequences for nwalign
    gapch = gapchars(j);
    subseq = subseq(subseq ~= gapch);
    if ~isempty(subseq)
        [subscore, bestalign, startat] = nwalign(subseq, params.seq{j}, 'alphabet', swalphas{j}, 'glocal', true);
    else
        subscore = 0;
        bestalign = [repmat('-', 1, length(params.seq{j})) ; blanks(length(params.seq{j})) ; params.seq{j}];
        startat = [1 ; 1];
    end
    aligndata{j}.score = subscore;
    aligndata{j}.alignment = bestalign;
    aligndata{j}.startat = startat;
end

% find the mapping to the reference sequence for each alphabet
for j = 1:nalphas
%    newalign.refseq(j).seq = params.seq{j};
    ngmap = findseqmap(aligndata{j}.alignment);
    % this is a little complicated because ngmap ignores the gaps
    subseq = alignment.data(params.idx, ranges(1, j):ranges(2, j));
    gapch = gapchars(j);
    map = zeros(length(subseq), 1);
    map(subseq ~= gapch) = ngmap;
    newalign.refseq(j).map = map;
end

% truncate, if needed
if params.truncate
    newalign.alphawidths = arrayfun(@(s) sum(s.map ~= 0), newalign.refseq);
    newdata = repmat(blanks(sum(newalign.alphawidths)), size(newalign.data, 1), 1);
    
    % check whether we need to remove any alphabets entirely
    mask = true(nalphas, 1);
    crt = 1;
    for j = 1:nalphas
        if newalign.alphawidths(j) == 0
            % this alphabet will be removed
            mask(j) = false;
        end
        
        subdata = newalign.data(:, ranges(1, j):ranges(2, j));
        subdata = subdata(:, newalign.refseq(j).map ~= 0);
        
        crt_end = crt + newalign.alphawidths(j) - 1;
        newdata(:, crt:crt_end) = subdata;
        
        newalign.refseq(j).map = newalign.refseq(j).map(newalign.refseq(j).map ~= 0);
        newalign.refseq(j).map = newalign.refseq(j).map(:);
        
        crt = crt_end + 1;
    end
    
    newalign.alphabets = newalign.alphabets(mask);
    newalign.alphawidths = newalign.alphawidths(mask);
    newalign.refseq = newalign.refseq(mask);    
    newalign.data = newdata;
end

% handle alternate position names, if provided
if ~isempty(params.posnames)
    % have cell array with one cell per alphabet
    if ~iscell(params.posnames) || ischar(params.posnames{1})
        params.posnames = {params.posnames};
    end
    for j = 1:nalphas
        if ~isempty(params.posnames{j})
            % map the position mappings to the names
            newalign.refseq(j).map = params.posnames{j}(newalign.refseq(j).map);
            newalign.refseq(j).map = newalign.refseq(j).map(:);
        end
    end
end

simcounts = cellfun(...
    @(align) sum(align.alignment(2, :) == '|' & isletter(align.alignment(1, :)) & isletter(align.alignment(3, :))), ...
    aligndata);
totals = cellfun(...
    @(align) sum(isletter(align.alignment(1, :)) & isletter(align.alignment(3, :))), ...
    aligndata);
alignfrac = sum(simcounts)/sum(totals);
if alignfrac < 1 && params.verbose
    warning([mfilename ':imperfect'], ['The matching between the reference sequence and the alignment is not perfect (identity = ' num2str(alignfrac, 3) ').']);
    for j = 1:nalphas
        disp([newalign.alphabets{j} ':']);
        disp(aligndata{j}.alignment);
    end
end

end