function [gapstats, changed] = statsincludegaps(stats, structure)
% STATSINCLUDEGAPS Include gap entries in the fields of the statistics
% structure.
%   gapstats = STATSINCLUDEGAPS(stats) includes entries for gaps in the
%   fields of the statistics structure, and updates its alphabets field as
%   done by alphaincludegaps.
%
%   gapfreq1 = STATSINCLUDEGAPS(freq1, structure) includes entries for gaps
%   in the frequency vector freq1, based on the structure given in
%   'structure' (which should have fields 'alphabets' and 'alphawidths').
%
%   gapcmat = STATSINCLUDEGAPS(cmat, structure) includes entries for gaps in
%   the covariance matrix cmat, based on the structure given in
%   'structure'.
%
%   Note that the gap entries cannot be recovered for the pairwise
%   frequency matrix freq2 without also having freq1, so the last form of
%   the function cannot be used with freq2.
%
%   [..., changed] = STATSINCLUDEGAPS(...) returns 'changed' equal to true
%   if a change needed to be made, false otherwise.
%
% See also: ALPHAINCLUDEGAPS, INCLUDEGAPS.

% Tiberiu Tesileanu (2014)

if ~statscheck(stats)
    if ~isnumeric(stats) || (~isvector(stats) && ~ismatrix(stats))
        error([mfilename ':badstats'], 'The first argument should be a statistics structure or a numeric vector or matrix.');
    end
    if nargin < 2
        error([mfilename ':toofew'], 'When the first argument is numeric, there should be a second argument indicating its structure.');
    end
    [gapstructure, changed] = alphaincludegaps(structure);
    
    if changed
        binmap = getbinmap(structure);
        gapbinmap = getbinmap(gapstructure);
        
        if isvector(stats)
            freq1 = stats(:);
            newfreq1 = zeros(gapbinmap{end}(end), 1);
            
            for i = 1:length(binmap)
                idxsnew = gapbinmap{i};
                idxsold = binmap{i};
                
                if length(idxsnew) == length(idxsold)
                    % nothing to change
                    newfreq1(idxsnew) = freq1(idxsold);
                else
                    newfreq1(idxsnew(2:end)) = freq1(idxsold);
                    newfreq1(idxsnew(1)) = 1 - sum(freq1(idxsold));
                end
            end
            
            gapstats = newfreq1;
        else
            cmat = stats;
            newcmat = zeros(gapbinmap{end}(end));
            
            for i = 1:length(binmap)
                idxs1new = gapbinmap{i};
                idxs1old = binmap{i};
                
                if length(idxs1new) == length(idxs1old)
                    idxs1match = idxs1new;
                else
                    idxs1match = idxs1new(2:end);
                end
                for j = 1:length(binmap)
                    idxs2new = gapbinmap{j};
                    idxs2old = binmap{j};
                    
                    if length(idxs2new) == length(idxs2old)
                        idxs2match = idxs2new;
                    else
                        idxs2match = idxs2new(2:end);
                    end
                    
                    block = cmat(idxs1old, idxs2old);
                    newcmat(idxs1match, idxs2match) = block;
                    if length(idxs1old) ~= length(idxs1new)
                        newcmat(idxs1new(1), idxs2match) = -sum(block, 1);
                    end
                    if length(idxs2old) ~= length(idxs2new)
                        newcmat(idxs1match, idxs2new(1)) = -sum(block, 2);
                        if length(idxs1old) ~= length(idxs1new)
                            newcmat(idxs1new(1), idxs2new(1)) = -sum(newcmat(idxs1new(1), idxs2new(2:end)));
                        end
                    end
                end
            end
            
            gapstats = newcmat;
        end
    else
        gapstats = stats;
    end
else
    if nargin > 1
        error([mfilename ':toomany'], 'Second argument is only used when the first argument is a numeric vector or a matrix.');
    end
    % make the alphabets gappy
    [gapstats, changed] = alphaincludegaps(stats);
    if changed    
        % add gap entries to freq1, cmat, and freq2 (if it exists)
        gapstats.freq1 = statsincludegaps(gapstats.freq1, stats);
        gapstats.cmat = statsincludegaps(gapstats.cmat, stats);
        if isfield(stats, 'freq2')
            gapstats.freq2 = gapstats.cmat + gapstats.freq1(:)*gapstats.freq1(:)';
        end
    else
        gapstats = stats;
    end
end

end