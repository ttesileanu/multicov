function mimat = getmi(alignment, varargin)
% GETMI Get matrix of pairwise mutual information.
%   mimat = GETMI(alignment) calculates the mutual information between
%   every two columns of the character or binary alignment, and returns it
%   as a matrix.
%
%   mimat = GETMI(stats) uses a statistics-like structure for the
%   calculation. Stats should have fields 'freq1', 'freq2' or 'cmat', and
%   'alphabets' and 'alphawidths' to describe their structure. The 'type'
%   field is not checked.

% Tiberiu Tesileanu (2014)

if alncheck(alignment) || bincheck(alignment)
    stats = getstats(alnincludegaps(alignment), 'freq2', true);
elseif isstruct(alignment) && all(isfield(alignment, {'freq1', 'alphabets', 'alphawidths'}))
    if ~any(isfield(alignment, {'freq2', 'cmat'}))
        error([mfilename ':nofreq2'], 'No freq2 or cmat fields in input; need at least one of them.');
    end
    if ~isfield(alignment, 'freq2') && isfield(alignment, 'cmat')
        alignment.freq2 = alignment.cmat + alignment.freq1(:)*alignment.freq1(:)';
    end
    [~, needsgap] = alphaincludegaps(alignment);
    if needsgap
        if ~isfield(alignment, 'cmat')
            alignment.cmat = alignment.freq2 - alignment.freq1(:)*alignment.freq1(:)';
        end
        alignment.freq1 = statsincludegaps(alignment.freq1, alignment);
        alignment.cmat = statsincludegaps(alignment.cmat, alignment);
        alignment.alphabets = alphaincludegaps(alignment.alphabets);
        alignment.freq2 = alignment.cmat + alignment.freq1(:)*alignment.freq1(:)';
    end
	stats = alignment;
else
    error([mfilename ':badarg'], 'The argument should be an alignment or a statistics structure.');
end

binmap = getbinmap(stats);
mimat = blockapply(stats.freq2, ...
    @(m, i, j) getmiblock(m, stats.freq1(binmap{i}), stats.freq1(binmap{j})), ...
    stats, 'indices', true);

end

function mi = getmiblock(m, freqi, freqj)
% GETMIBLOCK Get MI for a block.
%   mi = GETMIBLOCK(m, freqi, freqj) calculates the MI for the pairwise
%   frequency matrix m. For the calculation to be correct, freqi must equal
%   sum(m, 2) and freqj must equal sum(m, 1).

fprod = freqi(:)*freqj(:)';
mask = (m > eps & fprod > eps);
mi = sum(sum(m(mask) .* log(m(mask) ./ fprod(mask))));

end