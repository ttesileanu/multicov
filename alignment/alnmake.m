function align = alnmake(data, alphabet, annotations)
% ALNMAKE Make an alignment structure.
%   align = ALNMAKE creates a new, empty alignment. This is a structure
%   with fields alphabets, alphawidths, data, seqw, refseq, annotations,
%   and type.
%
%   Alphabets is a cell array of strings, identifying the kinds of data
%   contained in the alignment.
%
%   Alphawidths is a numeric array giving the widths of each subsequence
%   corresponding to each alphabet.
%
%   Data is the raw alignment data, presented as a character matrix.
%
%   Seqw is an array of sequence weights for the alignment.
%
%   Refseq is a structure array containing fields 'seqdb', 'seqid', and
%   'map', identifying a reference sequence and a mapping from alignment
%   positions to this sequence for each alphabet. Seqdb and seqid are
%   strings containing identification for the database and ID of the
%   reference sequence. Refseq(i).map is an array of size alphawidths(i)
%   such that its kth element gives the reference sequence position for the
%   kth column in the ith alphabet. The position can either be numeric or a
%   string (but this should be consistent for each alphabet). Refseq can
%   have an additional field 'seq' giving the actual reference sequence.
%
%   Annotations is a cell array of strings.
%
%   Type is set to 'character', to indicate a non-binary alignment.
%
%   By default, refseq.seqdb is set to 'generic', refseq.seqid is set to '',
%   and refseq.map(i) is set to 1:alphawidths(i). Also, the sequence weights
%   seqw are set to 1.
%
%   align = ALNMAKE(data, alphabet) creates an alignment with the given
%   alphabet, containing the data (given as a matrix).
%
%   ALNMAKE(data, alphabet, annotations) uses the given annotations.
%
%   See also ALNADD, ALNCHECK.

% Tiberiu Tesileanu (2012-2013)

% set the alignment identifier
align.type = 'character';

% handle the first form of the function
if nargin == 0
    align.alphabets = {};
    align.alphawidths = [];
    align.data = '';
    align.seqw = [];
    align.refseq = struct('seqdb', {}, 'seqid', {}, 'map', {});
    align.annotations = {};
    return;
end

% check parameters
if ~ismatrix(data) || ~ischar(data)
    error([mfilename ':baddata'], 'The first argument should be a character matrix.');
end
if ~ischar(alphabet)
    error([mfilename ':badalpha'], 'The second argument should be an alphabet.');
end

% handle the sequence weights
[N, w] = size(data);
if nargin < 3
    annotations = repmat({''}, N, 1);
end

align.alphabets = {alphabet};
align.alphawidths = w;
align.data = data;
align.seqw = ones(N, 1);
align.refseq = struct('seqdb', {'generic'}, 'seqid', {''}, 'map', {(1:w)'});
align.annotations = annotations;

end
