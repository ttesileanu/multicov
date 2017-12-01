function checked = showcontacts(contacts, diststruct, varargin)
% SHOWCONTACTS Make a figure showing the contact predictions compared to a
% structure.
%   checked = SHOWCONTACTS(contacts, diststruct) makes a figure displaying
%   the contact predictions in comparison to those contacts inferred from a
%   distance structure. It returns the output of checkcontacts.
%
%   Options:
%    'bkgcolors' {colcommon, colnostruct, colnopred}
%       Background colors to use. 'Colcommon' is the color for portions of
%       the molecule where both alignment and structural information
%       exists; 'colnostruct' is used for places where there is no
%       structural information; and 'colnopred' is used for places where
%       there is no alignment information. These values can be either RGB
%       triplets or Matlab character color codes.
%       (default: {'w', [0.7 0.7 1], [0.7 0.7 0.7]})
%    'contiguous' <b>
%       When one of the automatic 'range' options are used (i.e., there is
%       no 'forcerange'), if 'contiguous' is true, then the range that is
%       displayed contains a contiguous range in refseq. This only applies
%       for numeric refseqs.
%       (default: true)
%    'focus' <v/c>
%       Which alphabets to focus on. This can be a cell array with two
%       elements, each of which being a vector of alphabet indices.
%       {[a1, ..., am], [b1, ..., bn]} means use alphabets a1 to am on the
%       x-axis, and b1 to bn on the y-axis. Alternatively, this can be a
%       single vector v, in which case the same set of alphabets is used
%       for both axes.
%       (default: all alphabets are displayed)
%    'forcerange' <c>
%       Cell array choosing which positions to include for each alphabet.
%       The positions are given in terms of reference sequence indices.
%       (default: choose positions based on 'range' option)
%    'marker' <s>
%       The kind of marker used to show predicted contacts.
%       (default: '.')
%    'markercolors' {coltrue, colfalse, colunknown}
%       Colors to use for a correct contact prediction, for an incorrect
%       one, and for a case where the turth is not known. Each color can be
%       an RGB triplet or a Matlab character color code.
%       (default: {'g', 'r', [0.5 0.5 0.5]})
%    'markersize' <x>
%       The size of the markers used to show predicted contacts.
%       (default: 20)
%    'range' <s>
%       Which positions to show. This can be
%        'common':  Only the positions that are common between the distance
%                   structure and the contact predictions.
%        'contacts':Only the positions from the contacts structure.
%        'struct':  Only the positions from the distance structure.
%        'union':   All the positions from the distance structure and those
%                   from the contacts structure.
%       (default: 'union')
%    'revy' <b>
%       Whether the y-axis is reversed (larger positions going down).
%       (default: true)
%    'structcolor' <x>
%       Color to use for representing contacts in the distance structure.
%       This can be a Matlab character code or an RGB triplet.
%       (default: [0.4 0.4 0.4])
%    'threshold' <x>
%       Distance threshold to use for considering a contact (as in
%       checkcontacts).
%       (default: 6)
%
% See also: CHECKCONTACTS.

% Tiberiu Tesileanu (2013-2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('bkgcolors', {'w', [0.7 0.7 1], [0.7, 0.7, 0.7]}, ...
    @(c) iscell(c) && isvector(c) && length(c) == 3);
parser.addParamValue('contiguous', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('focus', [], @(v) isvector(v) && (isnumeric(v) || (iscell(v) && length(v) == 2)));
parser.addParamValue('forcerange', [], @(c) iscell(c) && isvector(c));
parser.addParamValue('marker', '.', @(c) ischar(c) && isscalar(c));
parser.addParamValue('markercolors', {'g', 'r', [0.5, 0.5, 0.5]}, ...
    @(c) iscell(c) && isvector(c) && length(c) == 3);
parser.addParamValue('markersize', 20, @(x) isscalar(x) && isnumeric(x));
parser.addParamValue('range', 'union', @(s) ismember(s, {'common', 'contacts', 'struct', 'union'}));
parser.addParamValue('revy', true, @(b) islogical(b) && isscalar(b));
parser.addParamValue('structcolor', [0.4 0.4 0.4], @(c) ischar(c) || (isnumeric(c) && length(c) == 3));
parser.addParamValue('threshold', 6, @(x) isscalar(x) && isnumeric(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

% check the contacts against the structure
if isempty(params.threshold)
    cc_opts = {};
else
    cc_opts = {'threshold', params.threshold};
end
checked = checkcontacts(contacts, diststruct, cc_opts{:});

% transform colors to RGB
params.bkgcolors = cellfun(@torgb, params.bkgcolors(:), 'uniform', false);
params.markercolors = cellfun(@torgb, params.markercolors(:), 'uniform', false);
params.structcolor = torgb(params.structcolor);

nalphas = length(contacts.refseq);
if nalphas ~= length(diststruct.refseq)
    error([mfilename ':badrefseq'], 'The lengths of the refseq structures should agree between the contacts and the structure.');
end

% handle the various forms of the 'focus' option
if isempty(params.focus)
    params.focus = 1:nalphas;
end
if ~iscell(params.focus)
    params.focus = {params.focus, params.focus};
end

if isempty(params.forcerange)
    % we need to decide on ranges based on the 'range' option
    switch params.range
        case 'common'
            params.forcerange = arrayfun(@(i) ...
                intersect(contacts.refseq(i).map, diststruct.refseq(i).map, 'stable'), ...
                1:nalphas, 'uniform', false);
        case 'contacts'
            params.forcerange = arrayfun(@(r) r.map, contacts.refseq, 'uniform', false);
        case 'struct'
            params.forcerange = arrayfun(@(r) r.map, diststruct.refseq, 'uniform', false);
        case 'union'
            params.forcerange = arrayfun(@(i) ...
                union(contacts.refseq(i).map, diststruct.refseq(i).map, 'stable'), ...
                1:nalphas, 'uniform', false);
        otherwise
            error([mfilename ':badrange'], 'Unrecognized option for ''range''.');
    end
    
    % make things contiguous if required
    if params.contiguous
        % only do this to numeric indices
        numeric_mask = cellfun(@(c) ~isempty(c) && isnumeric(c), params.forcerange);
        params.forcerange(numeric_mask) = cellfun(@(v) ...
            min(v(v>0)):max(v), params.forcerange(numeric_mask), 'uniform', false);
    end
end

% process the matrix we get from structure
struct_rendered = singlematrender(...
    diststruct.distmat < params.threshold, ...
    diststruct.validmask, ...
    diststruct.refseq, ...
    params.forcerange);

% process the matrix we get from contacts
contacts_rendered = singlematrender(...
    contacts.matrix, ...
    true(size(contacts.matrix)), ...
    contacts.refseq, ...
    params.forcerange);

% mark false positives
contacts_rendered((contacts_rendered > 0) & (struct_rendered == 0)) = 2;

% make a color map:
colormap([...
    params.bkgcolors{1} ; ...   % both structure and alignment data present, no structure contact
    params.bkgcolors{2} ; ...   % alignment data present, but no structure
    params.bkgcolors{3} ; ...   % structure data present with no contact, but no alignment data
    params.bkgcolors{2} + params.bkgcolors{3} - 1 ; ... % no structure data, no alignment data
    params.structcolor ; ...    % both structure and alignment data present, and structure contact
    params.structcolor ...      % structure data present with contact, but no alignment data
]);

% make the final matrix
matrix = zeros(size(struct_rendered));
matrix(struct_rendered < 0 & contacts_rendered < 0) = 4;
matrix(struct_rendered < 0 & contacts_rendered >= 0) = 2;
matrix(struct_rendered == 0 & contacts_rendered >= 0) = 1;
matrix(struct_rendered > 0 & contacts_rendered >= 0) = 5;
matrix(struct_rendered > 0 & contacts_rendered < 0) = 6;
matrix(struct_rendered == 0 & contacts_rendered < 0) = 3;

% find which part of this matrix we're supposed to display
alphawidths_res = cellfun(@length, params.forcerange);
alpharanges_res = getalpharanges(alphawidths_res);
xidxs = cell2mat(arrayfun(@(a) alpharanges_res(1, a):alpharanges_res(2, a), params.focus{1}(:)', 'uniform', false));
yidxs = cell2mat(arrayfun(@(a) alpharanges_res(1, a):alpharanges_res(2, a), params.focus{2}(:)', 'uniform', false));

submatrix = matrix(xidxs, yidxs);
subcontacts_rendered = contacts_rendered(xidxs, yidxs);
substruct_rendered = struct_rendered(xidxs, yidxs);

% display it!
image(submatrix);

was_hold = ishold;
hold on;

% now display the contacts
[true_contacts_i, true_contacts_j] = ind2sub(size(subcontacts_rendered), ...
    find(subcontacts_rendered > 0 & substruct_rendered > 0));
[false_contacts_i, false_contacts_j] = ind2sub(size(subcontacts_rendered), ...
    find(subcontacts_rendered > 1));
[unknown_contacts_i, unknown_contacts_j] = ind2sub(size(subcontacts_rendered), ...
    find(subcontacts_rendered > 0 & substruct_rendered < 0));

plot_opts = {params.marker, 'markersize', params.markersize};
plot(true_contacts_i, true_contacts_j, plot_opts{:}, 'color', params.markercolors{1});
plot(false_contacts_i, false_contacts_j, plot_opts{:}, 'color', params.markercolors{2});
plot(unknown_contacts_i, unknown_contacts_j, plot_opts{:}, 'color', params.markercolors{3});

% XXX draw some labels!

if params.revy
    set(gca, 'ydir', 'reverse');
else
    set(gca, 'ydir', 'normal');
end

if ~was_hold
    hold off;
end

end

function resmat = singlematrender(m, valid, refseq, ranges)
% SINGLEMATRENDER Render a contact matrix matching given reference sequence
% ranges.
%   resmat = SINGLEMATRENDER(m, valid, refseq, ranges) renders a contact
%   matrix based on given ranges of reference sequence positions. The
%   contact matrix 'm' and the matrix showing the validity of the contacts
%   'valid' use indices whose mapping to the reference sequences is given by
%   the 'refseq' structure. The resulting matrix 'resmat' uses the mapping
%   given by the cell array 'ranges'. The entries of the resulting matrix
%   follow the following code: they are equal to -1 if they don't correspond
%   to a valid entry in m; otherwise they are equal to the corresponding
%   entry in m.
%
%   Entries can be invalid if the reference sequence index (from 'ranges')
%   does not appear in 'refseq', or if valid is false at that position.

% invalid entries are coded as -1
m(~valid) = -1;

nalphas = length(refseq);
alphawidths_orig = arrayfun(@(s) length(s.map), refseq);
alpharanges_orig = getalpharanges(alphawidths_orig);
found_mask = [];
orig_idxs = [];
for a = 1:nalphas
    crtfound_mask = ismember(ranges{a}, refseq(a).map);
    if iscell(ranges{a})
        res2orig = cellfun(@(s) find(strcmp(refseq(a).map, s), 1), ranges{a}(crtfound_mask));
    else
        res2orig = arrayfun(@(i) find(refseq(a).map == i, 1), ranges{a}(crtfound_mask));
    end
    found_mask = [found_mask ; crtfound_mask(:)]; %#ok<AGROW>
    orig_idxs = [orig_idxs ; res2orig(:) + alpharanges_orig(1, a) - 1]; %#ok<AGROW>
end

alphawidths_res = cellfun(@length, ranges);
res_size = sum(alphawidths_res);
resmat = -1*ones(res_size, res_size);
found_mask = logical(found_mask);
resmat(found_mask, found_mask) = m(orig_idxs, orig_idxs);

end