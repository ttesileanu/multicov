function s = formattable(data, varargin)
% FORMATTABLE Format a table into a character array.
%   s = FORMATTABLE(data) arranges the data into a character array. The
%   data can be either a numeric, logical, or character array, or a cell
%   array. The elements of the cell array should be logical, numeric,
%   vectors, matrices, or strings. If data is a cell array of cell arrays,
%   then all of these should have the same length, as they will represent
%   each column of the table. Each column of the data will be displayed as
%   one formatted column in the final output.
%
%   Note that this can get quite slow, so it is not appropriate for using
%   with large datasets.
%
%   The following options control the output of this function:
%    'fixsize' <n/v>
%       Fix column sizes. If this is a single number, it applies to all
%       columns equally. If it is a vector, it should have one element per
%       column. A value of zero (per table or per column) means autodetect.
%       (default: 0)
%    'header' <c>
%       A cell array of strings to be used as a header. (default: none)
%    'hpos' <c/s>
%       Controls the horizontal positioning of elements in the columns.
%       This can be a single character (in which case the same choice
%       applies to all columns), or a string with a character for each
%       column. Choices:
%         'l': flushed left
%         'c': centered
%         'r': flushed right
%       (default: 'r')
%    'minsize' <x>
%       Set the minimum size of a column. (default: 3)
%    'maxsize' <x>
%       Set the maximum size of a column. When this is reached, the text in
%       each column is wrapped around. (default: 64)
%    'padding' <c/s>
%       Selects the character used for padding. This can be a single
%       character that applies for all columns, or a string containing one
%       character per column. (default: ' ')
%    'spacing' <n/v>
%       Number of spaces between columns. This can be a vector of length
%       (number of columns - 1) to have different spacings for different
%       columns. The entry v(i) corresponds to the space placed after the
%       ith column. (default: 2)
%    'vpos' <c>
%       Controls the horizontal positioning of elements in cases in which
%       a row of the table spans several lines of text. This can be
%         't': flushed to the top
%         'c': centered
%         'b': flushed to the bottom
%       (default: 'c')
%
%   See also: SPLITSTRING.

% Tiberiu Tesileanu (2013-2014)

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('fixsize', 0, @(x) isnumeric(x) && isvector(x));
parser.addParamValue('header', {}, @(c) iscell(c) && isvector(c));
parser.addParamValue('hpos', 'r', @(s) ischar(s) && isvector(s));
% XXX I don't know why minsize/maxsize and vpos can't be set per-column
parser.addParamValue('minsize', 3, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('maxsize', 64, @(x) isnumeric(x) && isscalar(x));
parser.addParamValue('padding', ' ', @(s) ischar(s) && isvector(s));
parser.addParamValue('spacing', 2, @(x) isnumeric(x) && isvector(x));
parser.addParamValue('vpos', 'c', @(c) ischar(c) && isscalar(c));

% parse
parser.parse(varargin{:});
params = parser.Results;

% extend 'fixsize', 'padding', 'hpos', and 'spacing' to one value per column
if length(params.fixsize) == 1
    params.fixsize = repmat(params.fixsize, 1, size(data, 2));
elseif length(params.fixsize) ~= size(data, 2)
    error([mfilename ':badfix'], 'Column size vector has size incompatible with data.');
end
if length(params.padding) == 1
    params.padding = repmat(params.padding, 1, size(data, 2));
elseif length(params.padding) ~= size(data, 2)
    error([mfilename ':badpad'], 'Padding vector has size incompatible with data.');
end
if length(params.hpos) == 1
    params.hpos = repmat(params.hpos, 1, size(data, 2));
elseif length(params.hpos) ~= size(data, 2)
    error([mfilename ':badhpos'], 'Horizontal position vector has size incompatible with data.');
end
if length(params.spacing) == 1
    params.spacing = repmat(params.spacing, 1, size(data, 2) - 1);
elseif length(params.spacing) ~= size(data, 2) - 1
    error([mfilename ':badsp'], 'Spacing vector has size incompatible with data.');
end
params.spacing = [params.spacing 0];

% the first thing to do: convert all the data to strings
% unfortunately this will be rather slow in general, since we want to allow
% any variety of data types

% handle non-cell-array data
if isnumeric(data) || islogical(data) || ischar(data)
    data = num2cell(data);
end
if ~iscell(data)
    error([mfilename ':baddata'], 'Unrecognized data format.');
end

% is this a cell array of cell arrays (columns)?
if isvector(data)
    % check that all columns are the same height
    if any(cellfun(@length, data) ~= length(data{1}))
        error([mfilename ':badcoldata'], 'When data is given as a cell array of columns, all columns should have the same height.');
    end
    % convert it to a matrix
    newdata = cell(length(data{1}), length(data));
    for j = 1:size(newdata, 2)
        crtcol = data{j};
        if ~iscell(crtcol)
            crtcol = num2cell(crtcol);
        end
        newdata(:, j) = crtcol(:);
    end
    data = newdata;
end

% go through the data and convert it to strings
for i = 1:size(data, 1)
    for j = 1:size(data, 2)
        if isnumeric(data{i, j})
            data{i, j} = num2str(data{i, j});
        elseif islogical(data{i, j})
            map = {'false', 'true'};
            data{i, j} = map(data{i, j} + 1);
        elseif ~ischar(data{i, j})
            error([mfilename ':noconv'], ['Can''t convert element at position (' int2str(i) ', ' int2str(j) ') to a string.']);
        end
    end
end

% add the header to the data, if present
if ~isempty(params.header)
    data = [params.header ; data];
end

% now find the size of each column
colsizes = arrayfun(@(j) max(cellfun(@(c) size(c, 2), data(:, j))), 1:size(data, 2));
% fix the sizes that are fixed by the options
mask = (params.fixsize > 0);
colsizes(mask) = params.fixsize(mask);
% implement min and max sizes, except where the size was fixed
colsizes(~mask) = min(max(colsizes(~mask), params.minsize), params.maxsize);

% calculate width of final table
width = sum(colsizes) + sum(params.spacing);

% resize cells so they fit in their column spaces
% also perform the horizontal positioning
for j = 1:size(data, 2)
    sz = colsizes(j);
    for i = 1:size(data, 1)
        % first get rid of useless space
        item = strtrim(data{i, j});
        if size(item, 2) > sz
            % need some reformatting; do it line-by-line
            newitem = '';
            for k = 1:size(item, 1)
                newitem = [newitem ; splitstring(item(k, :), sz)]; %#ok<*AGROW>
            end
            item = newitem;
        end
        % now, line by line, get rid of useless space, and perform the
        % positioning
        newitem = repmat(blanks(sz), size(item, 1), 1);
        for k = 1:size(item, 1)
            % trim
            l = strtrim(item(k, :));
            % pad
            cpad = params.padding(j);
            switch params.hpos(j)
                case 'l'
                    l = [l repmat(cpad, 1, sz - length(l))];
                case 'r'
                    l = [repmat(cpad, 1, sz - length(l)) l];
                case 'c'
                    npad = (sz - length(l)) / 2;
                    
                    l = [repmat(cpad, 1, floor(npad)) l repmat(cpad, 1, ceil(npad))];
                otherwise
                    error([mfilename ':badposchar'], ['Unrecognized horizontal positioning option, ' params.hpos(j) '.']);
            end
            newitem(k, :) = l;
        end
        data{i, j} = newitem;
    end
end

% get the row sizes
% note that each row in the table might occupy more than one line in the
% final output
rowsizes = arrayfun(@(i) max(cellfun(@(c) size(c, 1), data(i, :))), 1:size(data, 1));
rowsizes(rowsizes == 0) = 1;

% populate the table row by row
s = repmat(blanks(width), sum(rowsizes), 1);
crtrow = 1;
for i = 1:size(data, 1)
    rowsz = rowsizes(i);
    crtcol = 1;
    for j = 1:size(data, 2)
        item = data{i, j};
        % take care of the vertical positioning, if needed
        if rowsz > 1
            % get rid of empty lines of text in item
            if size(item, 1) > 1
                mask = arrayfun(@(k) ~isempty(strtrim(item(k, :))), 1:size(item, 1));
                item = item(mask, :);
            end
            % fill up with rows until we have the correct height
            if size(item, 1) < rowsz
                emptyrow = blanks(size(item, 2));
                switch params.vpos
                    case 't'
                        item = [item ; repmat(emptyrow, rowsz - size(item, 1), 1)];
                    case 'b'
                        item = [repmat(emptyrow, rowsz - size(item, 1), 1) ; item];
                    case 'c'
                        nrows = (rowsz - size(item, 1))/2;
                        item = [repmat(emptyrow, floor(nrows), 1) ; item ; repmat(emptyrow, ceil(nrows), 1)];
                    otherwise
                        error([mfilename ':badposchar'], ['Unrecognized vertical positioning option, ' params.vpos '.']);
                end
            end
        end
        
        s(crtrow + (0:size(item, 1)-1), crtcol + (0:size(item, 2)-1)) = item;
   
        crtcol = crtcol + size(item, 2) + params.spacing(j);
    end
    crtrow = crtrow + rowsz;
end

end