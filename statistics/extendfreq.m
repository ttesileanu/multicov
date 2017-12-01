function [freq1, freq2] = extendfreq(alignment, alphafreq)
% EXTENDFREQ Extend per-alphabet frequencies to alignment-wide values.
%   [freq1, freq2] = EXTENDFREQ(alignment, alphafreq) uses the per-alphabet
%   frequency information provided in the cell array alphafreq to generate
%   a frequency vector for the whole alignment. To do this, alphafreq{i} is
%   essentially repeated alignment.alphawidths(i) times. When the second
%   output argument is used, a pairwise frequency matrix is also generated.
%   Note that the calculation for freq2 is only performed when this value
%   is used by the calling code.
%
%   Instead of vectors of frequencies, each entry in alphafreq can be one
%   of 'uniform' or 'database'. If equal to 'uniform', the frequency is
%   taken to be constant over all amino acids *including* the gap (though
%   note that the mock gap from gapped alphabets is ignored). If equal to
%   'database', the frequency is obtained using getdbbkg.
%
%   Alphafreq can also be a single string instead of a cell array, equal to
%   either 'uniform' or 'database'. In this case, the same option is
%   applied to all alphabets in the alignment.
%
%   When the alignment has a single alphabet, alphafreq can be directly the
%   vector of frequency, without the enclosing cell array.
%
% See also: GETDBBKG.

% Tiberiu Tesileanu (2014)

% handle non-cell array input
if ~iscell(alphafreq)
    if ischar(alphafreq)
        % apply to all
        alphafreq = repmat({alphafreq}, length(alignment.alphabets), 1);
    else
        alphafreq = {alphafreq};
    end
end
if length(alphafreq) ~= length(alignment.alphabets)
    error([mfilename ':badlen'], 'Number of elements in alphafreq is incompatible with number of alphabets in alignment.');
end
if ~all(cellfun(...
        @(c) (isnumeric(c) && isvector(c)) || ...
             (ischar(c) && ismember(c, {'uniform', 'database'})), ...
        alphafreq))
    error([mfilename ':badfreq'], 'The second argument should be a cell array of numeric vectors, or the strings ''uniform'' or ''database''.');
end

freq1 = cell2mat(arrayfun(...
    @(i) extendsingle(alignment.alphabets{i}, alphafreq{i}, alignment.alphawidths(i)), ...
    (1:length(alignment.alphabets))', 'uniform', false));

if nargout >= 2
    freq2 = freq1*freq1';
    binmap = getbinmap(alignment);
    for i = 1:length(binmap)
        idxs = binmap{i};
        freq2(idxs, idxs) = diag(freq1(idxs));
    end
end

end

function freq = extendsingle(alphabet, mode, n)
% EXTENDSINGLE Extend the given frequency vector (or the one generated from
% the option if mode is 'uniform' or 'database') by repeating it n times.

if ischar(mode)
    if strcmp(mode, 'uniform')
        L = length(alphagetletters(alphabet, 'nogap'));
        if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
            denom = L;
        else
            denom = L + 1;
        end
        mode = ones(L, 1) / denom;
    elseif strcmp(mode, 'database')
        mode = getdbbkg(alphabet);
    end
end

freq = repmat(mode(:), n, 1);

end