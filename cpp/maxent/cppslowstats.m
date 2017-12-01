% CPPSLOWSTATS Internal function.
%
% Evaluate statistics given maximum entropy model parameters using a slow
% method that sums over all the states. (C++ code)
%   [freq1, freq2, Z] = CPPSLOWSTATS(alphabetletters, alphawidths, couplings, params)
%   evaluates the single-site ('freq1') pairwise ('freq2') frequencies, as
%   well as the partition function ('Z') for the maximum entropy model with
%   parameters given by 'couplings'. This is done by summing over all
%   possible states. Alphabetletters is a cell array giving the alphabets to
%   use -- note that this contains the actual alphabet letters (including
%   the gap), i.e., the results from alphagetletters. Alphawidths is the
%   number of columns used for each alphabet. Couplings is the matrix of
%   couplings and fields, gaps included. Finally, params is a structure that
%   may contain any number of the following fields:
%    'dispinterval': time interval in seconds between progress displays
%    'verbose':      whether to output progress information or not
%
% See also: ALPHAGETLETTERS, GETMAXENTENERGIES.

% Tiberiu Tesileanu (2013-2014)
