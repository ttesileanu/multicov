function colrgb = torgb(col0)
% COLRGB Internal function.
%
% Transform color to RGB triplet.
%   colrgb = TORGB(col0) returns an RGB triplet for the color col0, which
%   can be either a Matlab character color code, or an RGB triplet (in this
%   latter case, it is returned as-is).
%
%   This can also be used for a string of characters, in which case the
%   result is a matrix.

% Tiberiu Tesileanu (2014)

if ischar(col0)
    colrgb = zeros(length(col0), 3);
    for i = 1:length(col0)
        switch col0(i)
            case 'y'
                colrgb(i, :) = [1 1 0];
            case 'm'
                colrgb(i, :) = [1 0 1];
            case 'c'
                colrgb(i, :) = [0 1 1];
            case 'r'
                colrgb(i, :) = [1 0 0];
            case 'g'
                colrgb(i, :) = [0 1 0];
            case 'b'
                colrgb(i, :) = [0 0 1];
            case 'w'
                colrgb(i, :) = [1 1 1];
            case 'k'
                colrgb(i, :) = [0 0 0];
            otherwise
                error([mfilename ':badcol'], ['Unrecognized color ''' col0(i) '''.']);
        end
    end
else
    colrgb = col0;
end

end