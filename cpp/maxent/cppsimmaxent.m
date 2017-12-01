% CPPSIMMAXENT Internal function.
%
% Run MCMC simulation of maximum entropy model. (C++ code)
%   results = CPPSIMMAXENT(iniseq, alphabetletters, alphawidths, couplings,
%                          nsteps, params)
%   runs a Markov-chain Monte Carlo simulation for the maximum entropy
%   model with the given couplings, starting at the given initial sequence.
%   Alphabetletters is a cell array giving the alphabets to use -- note that
%   this contains the actual alphabet letters (including the gap), i.e., the
%   results from alphagetletters. Alphawidths is the number of columns used
%   for each alphabet. Couplings is the matrix of couplings and fields, gaps
%   included. Nsteps is the number of steps to make; set to 0 to run
%   forever (until CTRL+C is pressed). Note that CTRL+C is supported
%   gracefully (i.e., all the results obtained until the interruption are
%   saved and returned). Finally, params is a structure that may contain
%   any number of the following fields:
%    'burnin':       number of steps used as burnin (no sequences or
%                    energies saved)
%    'dispinterval': time interval in seconds between progress displays
%    'recstep':      step to use for recording sequences after the burnin
%    'verbose':      whether to output progress information or not
%
%   The output is also a structure, with fields:
%    'acceptances': average acceptance rates for the rec_step steps up to
%                   each stored sequence
%    'energies':    the energies of the stored sequences
%    'nsteps':      number of steps completed; this can be smaller than the
%                   number requested if CTRL+C is used
%    'sequences':   a matrix holding the sequences stored every rec-step
%                   steps after the burnin.
%
% See also: ALPHAGETLETTERS, GETMAXENTENERGIES.

% Tiberiu Tesileanu (2013-2014)