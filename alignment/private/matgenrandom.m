function mat = matgenrandom(N, freq, alphabet)
% MATGENRANDOM Internal function.
%
% Generate random alignment matrix.
%   mat = MATGENRANDOM(N, freq, alphabet) generates a random alignment
%   matrix of length N, using the given alphabet, and the given frequency
%   distribution.
%
%   See also: ALNGENRANDOM.

% Tiberiu Tesileanu (2012-2014)

freq = freq(:)';
letters = alphagetletters(alphabet);
lettersnogap = alphagetletters(alphabet, 'nogap');
nletts = length(lettersnogap);
width = length(freq) / nletts;
mat = repmat(blanks(width), N, 1);
localeps = 0.001; % don't worry about small  differences
for i = 1:width
    idxs = ((i - 1)*nletts + 1):(i*nletts);
    
    localsum = sum(freq(idxs));
    if localsum > 1
        if localsum - 1 > localeps
            warning([mfilename ':badfreq'], ['Frequencies for residue sum to >1! (' sum(freq(idxs)) ')']);
        end
        subfreq = [0 freq(idxs)/localsum];
    else
        subfreq = [1-localsum freq(idxs)];
    end
    mat(:, i) = gencolumn(N, subfreq, letters);
end

end

function col = gencolumn(N, freq, letters)
% GENCOLUMN Generate a column with letters occurring with given frequencies.
%   col = GENCOLUMN(N, freq, letters) generates a column with letters
%   sampled at the given frequencies.

counts = floor(N*freq);

N0 = sum(counts);

col = blanks(0);
for i = 1:length(letters)
    col = [col repmat(letters(i), 1, counts(i))]; %#ok<AGROW>
end

% the problem is, counts won't sum up to exactly N... we'll use randsample
% to sample the remaining elements
col = [col randsample(letters, N - N0, true, freq)];

% now randomize
col = col(randperm(N));

% make it a column vector
col = col';

end