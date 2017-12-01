function [ungapstats, changed] = statsexcludegaps(stats, structure)
% STATSEXCLUDEGAPS Exclude gap entries from the fields of the statistics
% structure.
%   gapstats = STATSEXCLUDEGAPS(stats) excludes entries for gaps from the
%   fields of the statistics structure, and updates its alphabets field as
%   done by alphaexcludegaps.
%
%   ungapfreq1 = STATSEXCLUDEGAPS(freq1, structure) excludes entries for gaps
%   from the frequency vector freq1, based on the structure given in
%   'structure' (which should have fields 'alphabets' and 'alphawidths').
%
%   ungapcmat = STATSEXCLUDEGAPS(cmat, structure) excludes entries for gaps
%   from  the covariance matrix cmat, based on the structure given in
%   'structure'.
%
%   This latter form can also be applied to remove the gap entries from
%   the pairwise frequency matrix.
%
%   [..., changed] = STATSEXCLUDEGAPS(...) returns 'changed' equal to true
%   if a change needed to be made, false otherwise.
%
% See also: ALPHAEXCLUDEGAPS, EXCLUDEGAPS, STATSINCLUDEGAPS.

% Tiberiu Tesileanu (2014)

if ~statscheck(stats)
    if ~isnumeric(stats) || (~isvector(stats) && ~ismatrix(stats))
        error([mfilename ':badstats'], 'The first argument should be a statistics structure or a numeric vector or matrix.');
    end
    if nargin < 2
        error([mfilename ':toofew'], 'When the first argument is numeric, there should be a second argument indicating its structure.');
    end
    [ungapstructure, changed] = alphaexcludegaps(structure);
    
    if changed
        binmap = getbinmap(structure);
        ungapbinmap = getbinmap(ungapstructure);
        
        if isvector(stats)
            freq1 = stats(:);
            newfreq1 = zeros(ungapbinmap{end}(end), 1);
            
            for i = 1:length(binmap)
                idxsnew = ungapbinmap{i};
                idxsold = binmap{i};
                
                if length(idxsnew) == length(idxsold)
                    % nothing to change
                    newfreq1(idxsnew) = freq1(idxsold);
                else
                    newfreq1(idxsnew) = freq1(idxsold(2:end));
                end
            end
            
            ungapstats = newfreq1;
        else
            cmat = stats;
            newcmat = zeros(ungapbinmap{end}(end));
            
            for i = 1:length(binmap)
                idxs1new = ungapbinmap{i};
                idxs1old = binmap{i};
                
                if length(idxs1new) == length(idxs1old)
                    idxs1match = idxs1old;
                else
                    idxs1match = idxs1old(2:end);
                end
                for j = 1:length(binmap)
                    idxs2new = ungapbinmap{j};
                    idxs2old = binmap{j};
                    
                    if length(idxs2new) == length(idxs2old)
                        idxs2match = idxs2old;
                    else
                        idxs2match = idxs2old(2:end);
                    end
                    
                    block = cmat(idxs1match, idxs2match);
                    newcmat(idxs1new, idxs2new) = block;
                end
            end
            
            ungapstats = newcmat;
        end
    else
        ungapstats = stats;
    end
else
    if nargin > 1
        error([mfilename ':toomany'], 'Second argument is only used when the first argument is a numeric vector or a matrix.');
    end
    % make the alphabets ungappy
    [ungapstats, changed] = alphaexcludegaps(stats);
    if changed    
        % add gap entries to freq1, cmat, and freq2 (if it exists)
        ungapstats.freq1 = statsexcludegaps(ungapstats.freq1, stats);
        ungapstats.cmat = statsexcludegaps(ungapstats.cmat, stats);
        if isfield(stats, 'freq2')
            ungapstats.freq2 = ungapstats.cmat + ungapstats.freq1(:)*ungapstats.freq1(:)';
        end
    else
        ungapstats = stats;
    end
end

end