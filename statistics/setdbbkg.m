function setdbbkg(alphabet, bkgfreq)
% SETDBBKG Set an alternative background frequency for the given alphabet.
%   SETDBBKG(alphabet, bkgfreq) replaces the background frequency for the
%   given alphabet. Bkgfreq should be a vector of numbers as long as the
%   number of letters in the alphabet, and these numbers should sum to 1.
%
%   SETDBBKG(alphabet, []) reverts to the default frequencies.
%
% See also: GETDBBKG.

% Tiberiu Tesileanu (2014)

global multicov_bkgfreq_database;

if isempty(multicov_bkgfreq_database)
    multicov_bkgfreq_database = containers.Map;
end

if ~isempty(bkgfreq)
    multicov_bkgfreq_database(alphabet) = bkgfreq(:);
else
    if multicov_bkgfreq_database.isKey(alphabet)
        multicov_bkgfreq_database.remove(alphabet);
    end
end

end