function bkgfreq = getdbbkg(alphabet, varargin)
% GETDBBKG Get background frequency for the given alphabet.
%   bkgfreq = GETDBBKG(alphabet) gets the background frequency of letters
%   in the alphabet, averaged over some large database. For proteins, these
%   are the values used for SCA. For other alphabets, the distribution that
%   is returned is currently uniform. See setdbbkg for changing the
%   distributions returned by this function.
%
%   Gapped alphabets return the same distribution as non-gapped alphabets
%   (setting the background frequency for gaps to exactly 0). This is true
%   even if the distribution of the non-gapped alphabet has been changed by
%   setdbbkg. It can, however, be changed by using setdbbkg on the gapped
%   alphabet.
%
% See also GETFREQ1, APPLYPOSW, PWFCTDKL, SETDBBKG.

% Tiberiu Tesileanu (2012-2014)

global multicov_bkgfreq_database;

% allow the user to override the defaults
if isempty(multicov_bkgfreq_database)
    multicov_bkgfreq_database = containers.Map;
else
    if multicov_bkgfreq_database.isKey(alphabet)
        bkgfreq = multicov_bkgfreq_database(alphabet);
        return;
    end
end

% gapped alphabets should have the same distribution as non-gapped ones,
% unless explicitly changed with setdbbkg
if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
    bkgfreq0 = getdbbkg(alphabet(4:end));
    bkgfreq = [0 bkgfreq0(:)'];
    bkgfreq = bkgfreq(:);
    return;
end

if strcmp(alphabet, 'protein')
    bkgfreq = [.073 .025 .050 .061 .042 .072 .023 .053 .064 .089 .023 .043 .052 .040 .052 .073 .056 .063 .013 .033];
else
    % for all others alphabets use a uniform distribution
    n = length(alphagetletters(alphabet, 'nogap'));
    bkgfreq = ones(1, n) / n;
end

bkgfreq = bkgfreq(:);

end