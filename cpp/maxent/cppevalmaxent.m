% CPPEVALMAXENT Internal function.
%
% Evaluate maximum entropy model energies for a set of sequences. (C++ code)
%   energies = CPPEVALMAXENT(sequences, alphabetletters, alphawidths, couplings)
%   evaluates the maximum entropy model energies for the given sequences.
%   Alphabetletters is a cell array giving the alphabets used for the
%   sequences -- note that this contains the actual alphabet letters
%   (including the gap), i.e., the results from alphagetletters.
%   Alphawidths is the number of columns used for each alphabet. Couplings
%   is the matrix of couplings and fields, gaps included.
%
% See also: ALPHAGETLETTERS, GETMAXENTENERGIES.

% Tiberiu Tesileanu (2013-2014)