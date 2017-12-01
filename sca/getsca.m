function stats = getsca(alignment, varargin)
% GETSCA Calculate SCA matrix.
%   stats = GETSCA(alignment) calculates the SCA matrix and returns it
%   together with other statistical information in the stats structure.
%   This works with either character or binary alignments. By default, it
%   uses the reduction method of SCA; see the options below for other
%   choices.
%
%   stats = GETSCA(stats0) can be used to calculate the SCA matrix from a
%   precalculated statistics structure.
%
%   Note that the positional weights are applied to all the fields of the
%   returned stats structure (freq1, cmat, and freq2, if it exists).
%
%   The return value of the function is a structure with fields
%    'alphabets'
%    'alphawidths'  -- duplicates the corresponding fields from the alignment
%    'cmat'         -- the SCA matrix
%    'cmatfull'     -- the SCA matrix before reduction; this is not present
%                      if no reduction function was specified, or if the
%                      'full' option is false
%    'freq1'        -- the single-site frequencies
%    'freq2'        -- the pairwise frequencies; this is not present if the
%                      'freq2' option is false
%    'type'         -- equal to 'stats'
%
%   Options:
%     'abs' <b>
%       Whether to take the absolute value of the covariance matrix.
%       (default: true)
%     'freq2' <b>
%       Whether to include the matrix of pairwise frequencies in the return
%       structure.
%       (default: false, to conserve space)
%     'full' <b>
%       Whether to also return the unreduced covariance matrix when
%       reduction is used.
%       (default: false, to conserve space)
%     'method' <s>
%       Choose a method for calculating the SCA matrix. Options are
%        'reduction':   the reduction method from SCA 4.0
%        'binary':      the binary approximation method from Halabi et al.
%        'projection':  the projection method from SCA 5.0
%       (default: 'reduction')
%     'pwfct' <f>
%       Specify the function that is used to calculate the positional
%       weights. Check applyposw for a description of the functions that
%       are allowed. Set this to an empty matrix to not perform any
%       positional weighting.
%       (default: @pwfctdkl, SCA's standard positional weighting function)
%     'redfct' <f>
%       Specify a reduction function that is applied to each residue block
%       of the covariance matrix and returns a number. This is commonly the
%       Frobenius or spectral norm. Set this to an empty matrix to not
%       perform the reduction. Note that this is ignored unless 'method' is
%       set to 'reduction'.
%       (default: @(m) norm(m, 'fro'))
%
%   See also GETSTATS, APPLYPOSW, PWFCTDKL, BLOCKAPPLY.

% Tiberiu Tesileanu (2012-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('abs', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('freq2', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('full', false, @(b) isscalar(b) && islogical(b));
parser.addParamValue('method', 'reduction', @(s) ismember(s, {'reduction', 'binary', 'projection'}));
parser.addParamValue('pwfct', @pwfctdkl, @(p) isempty(p) || isa(p, 'function_handle'));
parser.addParamValue('redfct', @(m) norm(m, 'fro'), @(p) isempty(p) || isa(p, 'function_handle'));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~alncheck(alignment) && ~bincheck(alignment) && ~statscheck(alignment)
    error([mfilename ':badaln'], 'First argument should either be a character or a binary alignment, or a statistics structure.');
end

if strcmp(params.method, 'reduction')
    % reduction method is easy
    if statscheck(alignment)
        stats = alignment;
    else
        stats = getstats(alignment, 'freq2', params.freq2);
    end
    stats = applyposw(stats, params.pwfct);
    if ~isempty(params.redfct)
        cmatfull = stats.cmat;
        if params.full
            stats.cmatfull = cmatfull;
        end
        stats.cmat = blockapply(cmatfull, params.redfct, stats);
    end
else
    % we do things slightly differently if we're given an alignment vs. a
    % statistics structure
    if ~statscheck(alignment)
        if alncheck(alignment)
            binalign = alntobin(alignment);
        else
            binalign = alignment;
        end
        freq1 = getfreq1(binalign);
        if ~isempty(params.pwfct)
            posw = params.pwfct(freq1, alignment);
        else
            posw = ones(size(freq1));
        end
        if strcmp(params.method, 'binary')
            % XXX what if pwamonut is not 1?!?!
            
            % now find which indices to keep -- the ones corresponding to
            % the consensus sequence
            considxs = getconsensus(binalign, 'indices', true);
            subposw = posw(considxs);
            subposw = subposw(:)';
            bindata = bsxfun(@times, full(binalign.data(:, considxs)), subposw);
            
            subfreq = freq1(considxs);
            stats.freq1 = subfreq(:).*subposw(:);
            freq2 = full(bindata'*diag(sparse(binalign.seqw(:)))*bindata/sum(binalign.seqw(:)));
            stats.cmat = freq2 - stats.freq1(:)*stats.freq1(:)';
            if params.freq2
                stats.freq2 = freq2;
            end
            stats.alphabets = alignment.alphabets;
            stats.alphawidths = alignment.alphawidths;
            stats.type = 'stats';
        elseif strcmp(params.method, 'projection')
            % XXX what if pwamonut is not 1?!?!
            
            binalign.data = bsxfun(@times, binalign.data, posw(:)');
            % need positionally-weighted frequencies
            binproj = binproject(binalign, getfreq1(binalign));
            
            stats = getstats(binproj);
        end
    else
        % we were given a statistics structure!
        % make sure we add or remove the freq2 field, as requested
        stats = getstats(alignment, 'freq2', params.freq2);
        
        if ~isempty(params.pwfct)
            posw = params.pwfct(stats.freq1, alignment);
        else
            posw = ones(size(stats.freq1));
        end
        
        if strcmp(params.method, 'binary')
            % find the indices to keep
            considxs = getconsensus(stats, 'indices', true);
            
            subfreq = stats.freq1(considxs);
            subposw = posw(considxs);
            subposw2 = subposw(:)*subposw(:)';
            
            stats.cmat = stats.cmat(considxs, considxs) .* subposw2;
            if params.freq2
                stats.freq2 = stats.freq2(considxs, considxs) .* subposw2;
            end
            stats.freq1 = subfreq(:) .* subposw;
        elseif strcmp(params.method, 'projection')
            stats = applyposw(stats, params.pwfct);
            projvec = stats.freq1;
            binmap = getbinmap(stats);
            for i = 1:length(binmap)
                nv = norm(projvec(binmap{i}));
                if nv > 0
                    projvec(binmap{i}) = projvec(binmap{i}) / nv;
                end
            end
            stats.freq1 = blockapply(stats.freq1(:) .* projvec(:), @sum, stats);
            stats.cmat = blockapply(stats.cmat .* (projvec(:)*projvec(:)'), @(m) sum(m(:)), stats);
            if params.freq2
                stats.freq2 = stats.cmat + stats.freq1(:)*stats.freq1(:)';
            end
        end
    end
end

% take the absolute value if requested
if params.abs
    stats.cmat = abs(stats.cmat);
end

end