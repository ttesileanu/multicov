function [dmodel, nfailed] = calculatedirectmodel(maxentparams, freq1, varargin)
% CALCULATEDIRECTMODEL Calculate the probability distributions of the
% 'direct model' employed in DCA.
%   dmodel = CALCULATEDIRECTMODEL(params, freq1) calculates the direct
%   model that is employed in DCA (see Weigt et al. 2009, Morcos et al. 2011).
%   The calculation uses the maximum entropy model couplings from the
%   structure 'params', and it also needs the frequency vector 'freq1'.
%   These parameters are required to 'include' their gaps, but CALCULATEDIRECTMODEL
%   calls the includegaps functions in case they don't. This works even if
%   one of them does include the gaps and the other doesn't.
%
%   The direct model is defined as follows. For each pair of sites (i, j), a
%   probability distribution P^{dir}_{ij}(A, B) is defined such that
%     P^{dir}_{ij}(A, B) = exp(ht_i(A) + ht_j(B) + J_{ij}(A, B)) ,
%   where J_{ij} is the corresponding block of couplings taken from 'params',
%   and ht_i(A) and ht_j(B) are calculated such that the marginals
%     P^{dir}_i(A) = \sum_j P^{dir}_{ij}(A, B) , and
%     P^{dir}_j(B) = \sum_i P^{dir}_{ij}(A, B)
%   match the values from 'freq1'. The values P^{dir}_{ij}(A, B) are then
%   returned in a matrix that has the same structure as the coupling matrix
%   (the one including gaps). Note that ht_i and ht_j depend on both sites
%   i and j.
%   The block-diagonal values (corresponding to P_{ii}(A, B)) are not
%   defined by this procedure, and are simply set to zero. As a consequence,
%   the block-diagonal couplings (i.e., the fields) are also not used in
%   this calculation.
%
%   Note that it is assumed that the coupling matrix is symmetric (and thus
%   in the output, P^{dir}_{ij}(A, B) = P^{dir}_{ji}(B, A)).
%
%   The returned structure has the following fields
%    'type'     -- set to 'directmodel'
%    'freq1'    -- equal to the single-site frequencies (including gaps)
%    'freq2'    -- equal to the values P^{dir}_{ij}(A, B) described above
%    'alphabets'
%    'alphawidths'
%    'refseq'
%               -- copied from params (except for promoting alphabets to
%                  gapped alphabets); describing the structure of freq1 and
%                  freq2
%
%   [dmodel, nfailed] = CALCULATEDIRECTMODEL(...) also returns the number
%   of pairs of sites at which the maximum number of iterations was reached
%   (i.e., the values for ht haven't converged to within the required
%   tolerance).
%
%   Options:
%    'dispinterval' <x>
%       The interval at which progress information is displayed (in seconds).
%       (default: 10)
%    'maxiter' <n>
%       Maximum number of iterations allowed per pair of sites to calculate
%       ht_i(A).
%       (default: 200)
%    'tol' <x>
%       Tolerance value to use for the iteration that calculates ht_i(A).
%       (default: 1e-4)
%    'verbose' <b>
%       Whether to display progress information and warnings (in case
%       'maxiter' is reached).
%       (default: true)
%
% See also: GETMAXENTPARAMS, INCLUDEGAPS.

% The code is adapted from the DCA code downloadable from Martin Weigt,
% http://dca.upmc.fr/DCA/Download.html
% (see also Weigt et al. 2009)

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('dispinterval', 10, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('maxiter', 200, @(n) isnumeric(n) && isscalar(n));
parser.addParamValue('tol', 1e-4, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('verbose', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

% make sure the parameters include the gaps
maxentparams = paramsincludegaps(maxentparams);
% make sure the frequencies include the gaps
binmap = getbinmap(maxentparams);
len = binmap{end}(end);
if length(freq1) ~= len
    freq1 = statsincludegaps(freq1, alphaexcludegaps(maxentparams));
end

freq1 = freq1(:);

dmodel.type = 'directmodel';
dmodel.freq1 = freq1;
dmodel.alphabets = maxentparams.alphabets;
dmodel.alphawidths = maxentparams.alphawidths;
dmodel.refseq = maxentparams.refseq;
dmodel.freq2 = zeros(size(maxentparams.couplings));

n = length(binmap);
nfailed = 0;
% XXX this would be much faster if we could vectorize it in some way
% for progress information: an approximation of the total number of
% elements to be calculated
hasdispd = false; % have we displayed any progress information?
timer = tic;
for i = 1:n
    if params.verbose
        dt = toc(timer);
        if dt >= params.dispinterval
            progress = (i - 1) * (2*n - i + 2) / n / (n + 1) * 100;
            disp(['DI calculation ' num2str(progress, 3) '% complete...']);
            drawnow('update');
            hasdispd = true;
            timer = tic;
        end
    end
    idxsi = binmap{i};
    Pi = freq1(idxsi);
    for j = (i + 1):n
        idxsj = binmap{j};
        Pj = freq1(idxsj);
        Jblock = maxentparams.couplings(idxsi, idxsj);
        [dmodel.freq2(idxsi, idxsj), failed] = getdirmodel(Jblock, Pi, Pj, params);
        nfailed = nfailed + failed;
    end
end
dmodel.freq2 = dmodel.freq2 + dmodel.freq2';
if hasdispd
    disp('DI calculation done.');
    drawnow('update');
end

end

function [Pdir, failed] = getdirmodel(Jblock, Pi, Pj, params)
% GETDIRMODEL Calculate direct probability model for a pair of sites.
%   [Pdir, failed] = GETDIRMODEL(Jblock, Pi, Pj, params) calculates the
%   direct probability model for the pair of sites, given the coupling
%   matrix Jblock and the marginals Pi and Pj. The function returns failed
%   = true if the maximum number of iterations (params.maxiter) is exceeded
%   without convergence to the given tolerance (params.tol).

% The code is adapted from the DCA code downloadable from Martin Weigt,
% http://dca.upmc.fr/DCA/Download.html
% (see also Weigt et al. 2009)

% Support for non-square blocks was added (no change was actually needed to
% the iteration code itself).

Ni = length(Pi);
Nj = length(Pj);
W = exp(Jblock);

% this code uses row vectors instead of column vectors
Pi = Pi(:)';
Pj = Pj(:)';

mui = ones(1, Ni) / Ni;
muj = ones(1, Nj) / Nj;

diff = 1.0;
iter = 0;
failed = false;
while diff > params.tol
	scrai = muj*W';
    scraj = mui*W;
    
	newi = Pi ./ scrai;
    newi = newi / sum(newi);
    
	newj = Pj ./ scraj;
    newj = newj / sum(newj);
    
	diff = max([max(abs(newi - mui)), max(abs(newj - muj))]);
    
	mui = newi;
    muj = newj;
    
    iter = iter + 1;
    if iter > params.maxiter
        if params.verbose
            warning([mfilename ':maxiter'], 'Maximum number of iterations reached for some pair of sites.');
        end
        failed = true;
        break;
    end
end

Pdir = W .* (mui'*muj);
% normalize the probabilities
Pdir = Pdir/sum(Pdir(:));

end