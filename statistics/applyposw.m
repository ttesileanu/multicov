function newstats = applyposw(stats, pwfct, pwamount)
% APPLYPOSW Apply positional weights to the statistics structure.
%   newstats = APPLYPOSW(stats, pwfct) applies the given positional
%   weight function to the statistics structure (as returned, for example,
%   by getstats). The function pwfct takes two arguments, the vector of
%   single-site statistics (as returned by getfreq1), and an alignment-like
%   structure, containing at least fields 'alphabets' and 'alphawidths'.
%   Note that some positional weighting functions may require additional
%   fields to be available in some cases. The positional weighting function
%   should return a vector w of positional weights of the same size as
%   freq1.
%
%   The positional weights are used as follows. Element (i, j) in the
%   covariance matrix is multiplied by w(i)*w(j). Element i in the
%   frequency vector is multiplied by w(i). The 'freq2' field, if it
%   exists, is updated to satisfy the relation between the covariance
%   matrix and freq1 and freq2.
%
%   APPLYPOSW(stats, pwfct, pwamount) performs an interpolation between
%   the positionally-weighted and the non-positionally-weighted statistics.
%   It simply returns stats if pwamount is zero.
%
% See also: GETSTATS.

% Tiberiu Tesileanu (2013-2014)

if nargin < 3
    pwamount = 1;
end

if pwamount > 0 && isa(pwfct, 'function_handle')
    posw = pwfct(stats.freq1, stats);
    newstats = stats;
    newstats.freq1 = posw(:).*stats.freq1(:);
    posw2 = posw(:)*posw(:)';
    newstats.cmat = posw2.*stats.cmat;
    
    if pwamount < 1 - eps
        newstats.freq1 = (1-pwamount)*stats.freq1 + pwamount*newstats.freq1;
        newstats.cmat = (1-pwamount)*stats.cmat + pwamount*newstats.cmat;
    end
    
    if isfield(stats, 'freq2')
        % XXX is this guaranteed to still have acceptable entries?
        % (i.e., >= 0, including for the implicit gap probabilities?)
        newstats.freq2 = newstats.cmat + newstats.freq1*newstats.freq1';
    end
else
    newstats = stats;
end

end