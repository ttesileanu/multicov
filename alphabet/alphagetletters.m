function letters = alphagetletters(alphabet, varargin)
% ALPHAGETLETTERS Get the letters in the alphabet.
%   letters = ALPHAGETLETTERS(alphabet) returns an ordered list of letters
%   in the given alphabet.
%
%   ALPHAGETLETTER(alphabet, 'nogap') does not include the gap in the
%   list. Otherwise, the gap will be the first letter in the list.
%
%   'Gapful' alphabets use a dummy character for a 'fake' gap, tricking
%   functions that usually ignore gaps to not ignore them; their name always
%   starts with the characters 'gap'. The fake gap character (usually a
%   blank, ' ') should not appear in the alignments for correct results.
%
%   list = ALPHAGETLETTERS('list') returns a cell array with all the
%   supported alphabets.

% Tiberiu Tesileanu (2012, 2014)

% handle the list of available alphabets
if nargin == 1 && strcmp(alphabet, 'list')
    % XXX this risks to fall out of sync with the actual code that supports
    % the alphabets
    letters = {'binary', 'dna', 'protein', 'rna'};
    letters = [letters arrayfun(@(i) ['multi' int2str(i)], 2:36, 'UniformOutput', false)];
    return;
end

% handle gapped alignments
if length(alphabet) > 3 && strcmp(alphabet(1:3), 'gap')
    % gapped version of alphabet
    letters = [' ' alphagetletters(alphabet(4:end))];
else
    recognized = true;
    % the standard alphabets
    switch alphabet
        case 'protein'
            letters = '-ACDEFGHIKLMNPQRSTVWY';
        case 'dna'
            letters = '-ACGT';
        case 'rna'
            letters = '-ACGU';
        case 'binary'
            letters = '01';
        otherwise
            recognized = false;
    end
    
    % the multiXX alphabets
    if ~recognized
        % there are a few more alphabets we can handle
        if length(alphabet) > 5 && strcmp(alphabet(1:5), 'multi')
            len = str2double(alphabet(6:end));
            if ~isnan(len) && isreal(len) && len == floor(len) && len <= 36 && len > 1
                recognized = true;
                ndigits = min(len, 10);                
                letters = arrayfun(@int2str, 0:(ndigits - 1));
                
                if len > 10
                    letters = [letters arrayfun(@(x) char('A' + x), 0:(len - 11))];
                end
            end
        end
    end
    
    if ~recognized
        error([mfilename ':badalpha'], ['Unrecognized alphabet ' alphabet '.']);
    end
end

% handle the 'nogap' option
if ~isempty(varargin)
    if strcmp(varargin{1}, 'nogap')
        letters = letters(2:end);
    else
        error([mfilename ':badopt'], 'The second argument, if provided, should be ''nogap''.');
    end
end

end