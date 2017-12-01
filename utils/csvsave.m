function csvsave(fname, table, varargin)
% CSVSAVE Save table to CSV file.
%   CSVSAVE(fname, table) saves the given table to a .csv file. The table
%   can be a numeric matrix, a 2d cell array, or a 1d cell array of
%   columns, each of which can be numeric or cell 1d arrays. Note that no
%   processing is performed on the data, so the file that is generated can
%   be invalid because of, e.g., delimiters or newlines present in fields.
%   Because of conversion and potential whitespace trimming, using CSVSAVE
%   after csvload or vice-versa will not do a perfect round-trip.
%
%   This function was designed to overcome the limitations of Matlab's
%   csvwrite, which forces the input to be a numeric matrix. However,
%   because CSVSAVE automatically detects the element types, it can be very
%   slow for large datasets.
%
%   Options:
%    'delimiter' <c>
%       What delimiter to use.
%       (default: ',')
%    'header' <c>
%       A 2d cell array representing a header to add to the table.
%       (default: {})
%
% See also: CSVLOAD.

% Tiberiu Tesileanu (2013-2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('delimiter', ',', @(c) ischar(c) && isscalar(c));
parser.addParamValue('header', {}, @(c) iscell(c) && isvector(c));

% parse
parser.parse(varargin{:});
params = parser.Results;

% identify the input type, and convert to 2d cell array
if ~iscell(table)
    if ~ismatrix(table)
        error([mfilename ':badndim'], 'The input table should be 2d or a 1d cell array of columns.');
    end
    table = num2cell(table);
else
    % if it's a 1d array, treat as array of columns
    if isvector(table)
        % make sure all columns have the same size
        vertical_size = length(table{1});
        newtable = cell(vertical_size, length(table));
        for i = 1:length(table)
            if length(table{i}) ~= vertical_size
                error([mfilename ':badshape'], 'All columns in the input data should have the same height.');
            end
            crtcol = table{i};
            if ~iscell(crtcol)
                crtcol = num2cell(crtcol);
            end
            newtable(:, i) = crtcol;
        end
        table = newtable;
    elseif ~ismatrix(table)
        error([mfilename ':badndim'], 'The input table should be 2d or a 1d cell array of columns.');
    end
end
% add header
if ~isempty(params.header) && size(params.header, 2) ~= size(table, 2)
    error([mfilename ':badhead'], 'Mismatch between number of columns in header and number of columns in table.');
end
table = [params.header ; table];
% make sure all elements are strings
for i = 1:size(table, 1)
    for j = 1:size(table, 2)
        if ~ischar(table{i, j})
            if iscell(table{i, j})
                error([mfilename ':badelem'], 'Table elements cannot be cell ararys.');
            end
            table{i, j} = num2str(table{i, j});
        end
    end
end

% write to file
f = fopen(fname, 'w');
for i = 1:size(table, 1)
    for j = 1:size(table, 2)
        if j > 1
            fprintf(f, params.delimiter);
        end
        fprintf(f, '%s', table{i, j});
    end
    fprintf(f, '\n');
end
fclose(f);

end