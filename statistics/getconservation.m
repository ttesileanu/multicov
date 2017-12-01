function cons = getconservation(alignment, varargin)
% GETCONSERVATION Get a conservation vector for the alignment.
%   cons = GETCONSERVATION(alignment) returns a vector of per-residue
%   conservation values for the alignment. By default, this is calculated
%   as the frequency of the consensus amino acid, excluding gaps. The
%   alignment can be in either character or binary form.
%
%   GETCONSERVATION(stats) finds the conservation using a statistics
%   structure instead. The only fields used are freq1, alphabets, and
%   alphawidths.
%
%   GETCONSERVATION(..., 'method', 'kl') defines conservation by the
%   KL divergence (relative entropy) between the empirical distribution of
%   amino acids at a site and a background distribution. The background
%   distribution for proteins is that used in work related to SCA, and is
%   obtained by averaging over a large library of proteins. For other
%   alphabets, the distribution is chosen to be uniform.
%
%   Note that when 'gaps' is true, the calculation of the KL divergence
%   requires a value for the background frequency of gaps. See the 'bkggap'
%   option for some choices related to this.
%
%   Available options are the following:
%     'bkggap' <str>/<x>/<vector>
%       Background gap frequency. This can be a character string, a number,
%       or a vector. If it is a string, it is one of:
%        'auto'
%           get the average gap frequency per alphabet from the current
%           alignment
%        'autoall'
%           get an average gap frequency for all alphabets in the current
%           alignment
%       If it is a number, this sets an overall gap frequency that is
%       constant for all positions. It can also be a vector with as many
%       elements as alphabets, giving a constant gap frequency for all the
%       residues within each alphabet; or a vector with as many elements as
%       there are positios in the alignment, giving a gap frequency for
%       each position.
%       (default: 'auto')
%     'method' <s>
%       This can be:
%        'kl':  calculate using the KL divergence
%        'max': calculate using the frequency of the consensus amino acid
%       (default: 'max')
%     'gaps' <b>
%       If this is false, gaps are ignored (i.e., they can't be consensus
%       amino acids for the 'max' method, and they are ignored for the 'kl'
%       method).
%       (default: false)

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('bkggap', 'auto', ...
    @(x) (ischar(x) && ismember(x, {'auto', 'autoall'})) || (isvector(x) && isnumeric(x)));
parser.addParamValue('method', 'max', @(s) ismember(s, {'max', 'kl'}));
parser.addParamValue('gaps', false, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

Npos = sum(alignment.alphawidths);

if isstruct(alignment) && all(isfield(alignment, {'freq1', 'alphabets', 'alphawidths'}))
    freq1 = alignment.freq1;
else
    if alncheck(alignment)
        binalign = alntobin(alignment);
    elseif bincheck(alignment)
        binalign = alignment;
    else
        error([mfilename ':badaln'], 'The first argument should be a character or binary alignment, or a statistics structure.');
    end
    freq1 = getfreq1(binalign);
end

% will be useful to know which columns are using gapped alphabets
gapped = false(Npos, 1);
ranges = getalpharanges(alignment.alphawidths);
for i = 1:length(alignment.alphabets)
    alphabet = alignment.alphabets{i};
    if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
        gapped(ranges(1, i):ranges(2, i)) = true;
    end
end

cons = zeros(Npos, 1);
alnmap = getbinmap(alignment);
if strcmp(params.method, 'max')
    for i = 1:Npos
        subfreq = freq1(alnmap{i});
        if params.gaps && ~gapped(i)
            subfreq = [1-sum(subfreq) ; subfreq(:)];
        elseif ~params.gaps && gapped(i)
            subfreq = subfreq(2:end);
        end
        cons(i) = max(subfreq);
    end
else
    if params.gaps
        % handle automatic background gap frequency determination
        if ischar(params.bkggap)
            gapfreq = zeros(length(alnmap), 1);
            gapfreq(~gapped) = cellfun(@(v) 1 - sum(freq1(v)), alnmap(~gapped));
            gapfreq(gapped) = cellfun(@(v) freq1(v(1)), alnmap(gapped));
            switch params.bkggap
                case 'auto'
                    % get gap frequency per each alphabet
                    params.bkggap = arrayfun(...
                        @(i) mean(gapfreq(ranges(1, i):ranges(2, i))), ...
                        1:length(alignment.alphabets));
                case 'autoall'
                    % get one global gap frequency
                    params.bkggap = mean(gapfreq);
                otherwise
                    error([mfilename ':badgap'], ['Unrecognized gap policy ''' params.bkggap '''.']);
            end
        end
        if isscalar(params.bkggap)
            % if we have a single number as global gap frequency:
            % extend it to a vector as wide as the alignment
            params.bkggap = repmat(params.bkggap, Npos, 1);
        elseif isvector(params.bkggap) && length(params.bkggap) == length(alignment.alphabets)
            % if we have one gap frequency per alphabet:
            % extend these to get a vector as wide as the alignment
            newbkggap = [];
            for i = 1:length(alignment.alphabets)
                crt = params.bkggap(i);
                newbkggap = [newbkggap ; repmat(crt, alignment.alphawidths(i), 1)]; %#ok<AGROW>
            end
            params.bkggap = newbkggap;
        elseif ~isvector(params.bkggap) || length(params.bkggap) ~= Npos
            error([mfilename ':badgap2'], ['The ''bkggap'' parameter should either be a string, a scalar, ' ...
                'or a vector as large as the number of alphabets or as wide as the alignment.']);
        end
    end
    
    bkgfreq1 = extendfreq(alignment, 'database');
    
    for i = 1:Npos
        block = freq1(alnmap{i});
        f = block(:);
        bkgf = flatten(bkgfreq1(alnmap{i}));
        if gapped(i)
            f = f(2:end);
            bkgf = bkgf(2:end);
        end

        if params.gaps
            % add a gap
            gapf = 1 - sum(f);
            f = [gapf ; f]; %#ok<AGROW>
            bkggapfreq = params.bkggap(i);
            bkgf = [bkggapfreq ; (1 - bkggapfreq)*bkgf];
        else
            s = sum(f);
            if s > 0
                f = f / s;
            end
        end
        
        mask = (f > eps);
        cons(i) = sum(f(mask) .* log(f(mask) ./ bkgf(mask)));
    end
end

cons = cons(:);

end