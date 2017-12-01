function res = pwfctdkl(freq, pattern)
% PWFCTDKL Calculate the standard SCA positional weights.
%   res = PWFCTDKL(freq, pattern) calculates the standard SCA positional
%   weights given the single-site statistics freq and the alignment
%   structure implied by pattern (which must contain fields alphabets and
%   alphawidths).
%
%   The weights are calculated as the KL divergence (relative entropy)
%   between freq and a set of background frequencies as used in SCA (which
%   were calculated by averaging over all known proteins at some time,
%   probably in 1999).
%
%   Note that whenever freq is exactly 0 or exactly 1, the relative entropy
%   diverges (since the background frequencies never reach these values).
%   To avoid problems, the values of freq that are <= 0 or >= 1 are
%   replaced by the most extreme allowed values present in the vector
%   (i.e., if freq(i) = 0, it is replaced by min(freq(freq > 0)), and if
%   freq(i) = 1, it is replaced by max(freq(freq < 1))). This function
%
% See also: APPLYPOSW.

% Tiberiu Tesileanu (2012, 2014)

if ~isvector(freq)
    error([mfilename ':badshape'], 'The first input should be a vector.');
end
if ~isstruct(pattern) || ~all(isfield(pattern, {'alphabets', 'alphawidths'}))
    error([mfilename ':badpatt'], 'The second input should be an alignment-like stucture.');
end

bkgfreq = extendfreq(pattern, 'database');

if length(freq) ~= length(bkgfreq)
    error([mfilename ':badlen'], 'The frequency vector doesn''t match the alignment structure.');
end

% have the result in the same shape as the input frequency vector
if size(freq, 1) == 1
    bkgfreq = bkgfreq';
end

% if all the frequencies are invalid, return zeros
if all(freq <= 0 | freq >= 1)
    res = zeros(size(freq));
else
    % fix the invalid frequencies
    mask = (freq <= 0);
    freq(mask) = min(freq(~mask));
    
    mask = (freq >= 1);
    freq(mask) = max(freq(~mask));
    
    % do the calculation
    res = log(freq.*(1 - bkgfreq) ./ ((1 - freq).*bkgfreq));
end

end