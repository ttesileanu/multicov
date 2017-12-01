function data = csvload(fname, varargin)
% CSVLOAD Load data from comma-separated value files.
%   s = CSVLOAD(fname) parses the .csv file, and returns the data. The data
%   is interpreted as a set of columns, which motivates the structure of
%   the output argument:
%     'header':   a cell array containing header rows; this may be empty if
%                 'header' is set to false, or if the automatic detection
%                 didn't find any header.
%     'columns':  the parsed data, organized as a cell array of columns.
%     'data':     the parsed data, organized as a structure with fieldnames
%                 that are related to column header names; this is only
%                 present if a header is found; columns whose header names
%                 cannot be transformed to field names get names of the
%                 form'col???'.
%     'coltypes': this is only present if 'process' is true, and it is a
%                 string identifying the column types: 'n' for numeric, 'v'
%                 for vectors, 's' for strings, 'l' for logical, 'x' for
%                 empty.
% 
%   If 'process' is false, all the data is returned in the form of strings.
%   Otherwise, each column in the data has a type given by coltypes: 'n'
%   columns are numeric arrays, 'v' columns are cell arrays of vectors, 'l'
%   columns are logical arrays, and 's' and 'x' columns are cell arrays of
%   strings.
%
%   This function was created to overcome the limitations of Matlab's
%   native .csv functions. The csvread and dlmread functions only work with
%   numeric data, while importdata behaves somewhat erratically (or at
%   least in an undocumented way) for files that mix numeric and
%   non-numeric data.
%
%   The automatic conversion facilities of this function are relatively
%   time-consuming, and can also use more memory than necessary, since the
%   whole file is first stored in memory in the form of a cell array of
%   strings. This function is thus not well-suited for very large
%   datasets, and even moderate datasets can be slow to load when
%   autodetection is required.
%
%   Options:
%     'comment' <c>
%       Character(s) representing a comment line. Set to empty to have no
%       comments. If there are several characters, each of them can be used
%       for comments.
%       (default: '#')
%     'delimiter' <s>
%       Delimiter to use.
%       (default: ',')
%     'header' <b>
%       If true, the top line is considered a header, and its entries are
%       used as column names. If false, there is no header, and the return
%       structure does not contain the 'data' field.
%       (default: autodetect)
%     'missingisnan' <b>
%       Whether to treat a missing number as NaN.
%       (default: true)
%     'namefct' <h>
%       Handle to function used to modify column names. After this is
%       applied, processing is still performed to make sure that the
%       returned name can be used as a fieldname in a structure. Set to
%       an empty array to perform no processing.
%       (default: @lower)
%     'process' <b>
%       Whether to process the data, or just return the it as raw cell
%       arrays of strings.
%       (default: true)
%     'trimws' <b>
%       Whether to trim whitespace at the edges of each cell.
%       (default: true)
%
% See also: CSVSAVE.

% Tiberiu Tesileanu (2013-2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('comment', '#', @(s) ischar(s) && isvector(s));
parser.addParamValue('delimiter', ',', @(c) ischar(c) && isscalar(c));
parser.addParamValue('header', [], @(b) isempty(b) || (isscalar(b) && islogical(b)));
parser.addParamValue('missingisnan', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('namefct', @lower, @(f) isempty(f) || isa(f, 'function_handle'));
parser.addParamValue('process', true, @(b) isscalar(b) && islogical(b));
parser.addParamValue('trimws', true, @(b) isscalar(b) && islogical(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if ~ischar(params.delimiter) || length(params.delimiter) > 1
    error([mfilename ':baddel'], 'Delimiter should be a single character.');
end

% first read everything into a cell table
f = fopen(fname);

rawdata = {};
i = 1;
while true
    line = fgetl(f);
    if ~ischar(line)
        % we're done -- fgetl returns -1 in case of EOF
        break;
    end
    
    % identify comments
    if ~isempty(params.comment)
        if params.trimws
            tmp = strtrim(line);
        else
            tmp = line;
        end
        if ismember(tmp(1), params.comment)
            % comment line, skip
            continue;
        end
    end
    
    % split the file up at the delimiters
    delpos = [0 find(line == params.delimiter) length(line)+1];
    items = arrayfun(@(i) line((delpos(i) + 1):(delpos(i+1)-1)), 1:(length(delpos)-1), 'uniform', false);
    if params.trimws
        items = cellfun(@strtrim, items, 'uniform', false);
    end
    
    % store the line
    rawdata(i, :) = items; %#ok<AGROW>
    i = i +1;
end

fclose(f);

% process this further only if we're asked for it, but there's nothing to do
% if the data part is empty
if ~params.process || isempty(rawdata) || (~isempty(params.header) && params.header && size(rawdata, 1) < 2)
    if isempty(params.header)
        params.header = false;
    end
    headerlines = double(params.header);
    data.header = rawdata(1:headerlines, :);
    converted = rawdata((headerlines+1):end, :);
else
    % try to find column types and autodetect header (if required)
    if ~isempty(params.header)
        headerlines = double(params.header);
        data.header = rawdata(1:headerlines, :);
        rawdata = rawdata((headerlines+1):end, :);
    else
        headerguesses = false(1, size(rawdata, 2));
    end
    coltypes = blanks(size(rawdata, 2));
    converted = rawdata;
    % get types for each element in the data file
    [types, values] = cellfun(@guesstype, rawdata, 'uniform', false);
    types = cell2mat(types);
    for i = 1:size(rawdata, 2)
        % get a list of all the types used in the column
        utypes = unique(types(:, i));
        % some conversions might be possible:
        % if there are v values around, we might as well turn all n's into v's
        if ismember('v', utypes) && ismember('n', utypes)
            mask = (types(:, i) == 'n');
            types(mask, i) = 'v';
            utypes = setdiff(utypes, 'n');
        end
        
        % also
        % x values can be converted to v, s, or n (if missignisnan is true)
        if length(utypes) > 1 && ismember('x', utypes)
            % if there are vectors around, or if missingisnan is false,
            % convert to them
            mask = (types(:, i) == 'x');
            if ismember('v', utypes) || ~params.missingisnan
                types(mask, i) = 'v';
                values{mask, i} = [];
            elseif ismember('n', utypes) && params.missingisnan
                types(mask, i) = 'n';
                values{mask, i} = nan;
            end
            % update the utypes
            utypes = setdiff(utypes, 'x');
        end

        if length(utypes) > 1 && isempty(params.header)
            % if we could get a single type provided there's a header, then
            % assume we have a header
            utypes0 = unique(types(2:end, i));
            if length(utypes0) == 1
                headerguesses(i) = true;
                utypes = utypes0;
            end
        end
        
        % if we still have more than one type in this column...
        if length(utypes) > 1
            coltypes(i) = 's';
        else
            coltypes(i) = utypes;
        end
        if coltypes(i) ~= 's' && coltypes(i) ~= 'x'
            converted(:, i) = values(:, i);
        end
    end
    if isempty(params.header)
        params.header = any(headerguesses);
        headerlines = double(params.header);
        data.header = rawdata(1:headerlines, :);
        converted = converted((headerlines+1):end, :);
    end
    
    data.coltypes = coltypes;
end

% present the data in column form
for j = 1:size(converted, 2)
    col = converted(:, j);
    if isfield(data, 'coltypes') && (data.coltypes(j) == 'n' || data.coltypes(j) == 'l')
        col = cell2mat(col);
    end
    data.columns{j} = col;
end
if ~isempty(data.header)
    % use header row to figure out names for columns
    names = data.header(end, :);
    % apply the user-specified transform function
    if ~isempty(params.namefct)
        names = cellfun(params.namefct, names, 'uniform', false);
    end
    % can't have names with no letters
    need_names = find(cellfun(@(s) all(~isletter(s)), names));
    for i = 1:length(need_names)
        names{need_names(i)} = ['col' int2str(need_names(i))];
    end
    % turn the names into things that are usable as field names
    names = cellfun(@tofieldname, names, 'uniform', false);
    % find duplicate names, and remove the ambiguity
    [unames, ~, j] = unique(names);
    if length(unames) < length(names)
        for k = 1:length(unames)
            m = sum(j == k);
            if m > 1
                % disambiguate
                places = find(j == k);
                for l = 1:length(places)
                    names{l} = [names{l} '_' int2str(l)];
                end
            end
        end
    end
    % finally, create the output structure
    for i = 1:length(names)
        data.data.(names{i}) = data.columns{i};
    end
end

end