function [alignment1, idx, oldidx] = seqsearch(alignment, varargin)
% SEQSEARCH Search for a particular sequence in the alignment.
%   alignment = SEQSEARCH(alignment, ...) searches for a particular
%   sequence in the alignment, and returns an alignment in which this
%   sequence is swapped to the first position. The sequence can be
%   identified using its index in the alignment (see 'idx' option), its
%   annotation (see 'annot' option), or its specific sequence (see 'seq'
%   option). Note that the sequence can be specified together with either
%   the index or the annotation, in which case the match between the given
%   sequence and the one found in the alignment is checked.
%
%   Note that sequence matching currently only works for 'protein', 'dna',
%   and 'rna' alphabets. The gapped versions ('gapprotein', 'gapdna',
%   'gaprna') behave identical to their ungapped counterparts.
%
%   [~, idx] = SEQSEARCH(alignment, 'swaptofirst', false, ...) keeps the
%   sequence where it is, and returns its index.
%
%   When 'swaptofirst' is true, the second return argument is always 1 (the
%   new position of the reference sequence), and there is a third return
%   argument that gives the position where the reference sequence was found
%   in the original alignment.
%
%   Options:
%    'annot' <s>
%       A regular expression used for matching the annotations in the
%       alignment in search for the sequence.
%    'idx' <i/v>
%       Identify the sequence by its position in the alignment. This can be
%       a vector of indices, in which case the actual sequence will be
%       identified using annotations or sequence, if provided.
%    'seq' <s/c>
%       This is either a string or a cell array of strings giving the
%       sequence that we are looking for. When this is a cell array, there
%       should be one string for each of the alphabets in the alignment.
%    'swaptofirst' <b>
%       Whether to swap the sequence into the first place in the returned
%       alignment (true), or to keep it where it is and return its index
%       (false).
%       (default: true)
%    'verbose' <b>
%       Whether to report ambiguous results to the user.
%       (default: true)

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('annot', [], @(s) ischar(s) && isvector(s));
parser.addParamValue('idx', [], @(i) isnumeric(i) && isvector(i) && isreal(i));
parser.addParamValue('seq', [], @(s) (ischar(s) && isvector(s)) || (iscell(s) && all(cellfun(@(ss) ischar(ss) && isvector(ss), s))));
parser.addParamValue('swaptofirst', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('verbose', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% check that the alignment is an alignment
if ~alncheck(alignment)
    error([mfilename ':badarg'], 'First argument should be an alignment structure.');
end

% what kind of search are we doing?
if ~isempty(params.idx)
    % we're actually allowing the index to be a range of indices
    idx = params.idx;
else
    idx = 1:size(alignment.data, 1);
end

if ~isempty(params.annot)
    % restrict to the sequences with the right annotations, if they are exist
    matches = (cellfun(@length, regexp(alignment.annotations, params.annot)) > 0);
    idx = intersect(idx, find(matches));
end

if ~isempty(params.seq)
    seq = params.seq;
    if ~iscell(seq)
        seq = {seq};
    end
    % finish up by using a sequence search
    if length(idx) > 1
        % ...but only if we have something to search through
        nalphas = length(alignment.alphabets);

        % first transform alphabets to a form that swalign/nwalign like
        swalphas = cell(1, nalphas);
        for i = 1:nalphas
            alphabet = alignment.alphabets{i};
            if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
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
        
        % next calculate the scores for each sequence
        scores = -inf(length(idx), 1);
        for i = 1:length(idx)
            ai = idx(i);
            score = 0;
            crtpos = 1;
            for j = 1:nalphas
                % get the subsequence for this alphabet
                crtpos_end = crtpos + alignment.alphawidths(j) - 1;
                subseq = alignment.data(ai, crtpos:crtpos_end);
                crtpos = crtpos_end;
                
                % get rid of gaps for the comparison
                subseq = subseq(subseq ~= '-');
                if ~isempty(subseq)
                    % nwalign requires non-empty sequences
                    score = score + nwalign(subseq, seq{j}, 'alphabet', swalphas{j}, 'glocal', true);
                end
            end
            scores(i) = score;
        end
        
        % finally, find the largest score(s)
        [~, subidx] = max(scores);
        idx = idx(subidx);
    end
end

% give out warnings if we need to
if params.verbose
    if isempty(idx)
        warning([mfilename ':notfound'], 'No sequence found.');
    else
        postext = 'position';
        if length(idx) > 1
            msg = 'Several matches found.';
            if params.swaptofirst
                msg = [msg ' Using the first one.'];
            end
            warning([mfilename ':manyfound'], msg);
            postext = 'positions';
        end
        disp(['Sequence found at ' postext ': ' int2str(idx) '.']);
    end
end

oldidx = idx;
alignment1 = alignment;

% swap the found sequence to first position, if asked to
if params.swaptofirst && ~isempty(oldidx)
    idx = 1;
    oldidx = oldidx(1);
    alignment1.data([idx oldidx], :) = alignment1.data([oldidx idx], :);
    alignment1.seqw([idx oldidx]) = alignment1.seqw([oldidx idx]);
    alignment1.annotations([idx oldidx]) = alignment1.annotations([oldidx idx]);
end

end