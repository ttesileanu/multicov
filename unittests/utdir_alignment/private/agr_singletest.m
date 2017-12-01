function testres = agr_singletest(expfreq, nreps, varargin)
% AGR_SINGLETEST Perform a single test for alngenrandom.
%   testres = AGR_SINGLETEST(expfreq, nreps, args...) generates nreps
%   alignments using alngenrandom with the arguments provided by args (the
%   arguments passed to AGR_SINGLETEST from the third position onwards),
%   and check that their distribution is consistent with the given expected
%   frequencies, expfreq. The function returns true if all the tests pass,
%   false otherwise.
%
%   Note that this function checks whether the generated distributions are
%   reasonably close to the expected ones; it is entirely possible that
%   this will fail by chance even if alngenrandom works properly. The
%   function also checks that alngenrandom doesn't always give the same
%   answers, which again it could do completely by chance (though very
%   rarely). For these reasons, care should be taken when using this
%   function.

% Tiberiu Tesileanu (2014)

% threshold for p value
pthresh = 0.05;

nfails = 0;
lastalign = [];
N = varargin{1};
for i = 1:nreps
    rndalign = alngenrandom(varargin{:});
    rndfreq = getfreq1(rndalign);
    
    failed = false;
    
    % we apply a chi^2 test
    % we need to take special care of cases in which the expected frequency
    % is zero
    zeromask = (expfreq < eps);
    if any(rndfreq(zeromask) > eps)
        failed = true;
    else
        chi2 = N*sum((rndfreq(~zeromask) - expfreq(~zeromask)).^2 ./ expfreq(~zeromask));
        p = 1 - chi2cdf(chi2, sum(~zeromask));
        if p < pthresh
            failed = true;
        end
    end

    % the generated alignments shouldn't all be identical
    if ~isempty(lastalign) && alncompare(lastalign, rndalign)
        failed = true;
    end
    
    if failed
        nfails = nfails + 1;
    end
    
    lastalign = rndalign;
end

expfails = pthresh*nreps;
stdfails = sqrt(expfails);
testres = (nfails <= (expfails + stdfails));

end