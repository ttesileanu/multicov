function checked = checkcontacts(contacts, diststruct, varargin)
% CHECKCONTACTS Check the accuracy of contact predictions.
%   checked = CHECKCONTACTS(contacts_pred, diststruct) returns an updated
%   version of the 'contacts_pred' structure that contains information about
%   the match between the contact predictions and the known contacts, as
%   inferred from the 'diststruct' structure. Contacts_pred should have the
%   format returned by getcontacts.
%
%   The returned structure is a copy of 'contacts_pred', with two extra
%   fields, 'check' and 'check_mapped'. 'Check' is a structure with fields
%    'correctmask'
%       A binary vector in which the ith entry shows whether the ith pair
%       of residues in contacts_pred.rawpairs matches a contact in
%       diststruct. See the 'threshold' option for the definition of a
%       contact.
%    'knownmask'
%       A binary vector in which the ith entry shows whether the distance
%       structure contains information for the ith pair of residues in
%       contacts_pred.rawpairs.
%    'precision'
%       A vector in which the ith entry shows the precision of the first i
%       contact predictions. Precision is defined as the fraction of
%       predicted contacts that are true contacts.
%    'tprate'
%       A vector in which the ith entry is the true positive rate for the
%       first i contact predictions. The true positive rate is the ratio
%       between the number of predicted contacts that are true contacts and
%       the total number of true contacts.
%    'fprate'
%       A vector in which the ith entry is the false positive rate for the
%       first i contact predictions. The false positive rate is the ratio
%       between the number of predicted contacts that are not true contacts
%       and the total number of pairs that are not contacts in 'diststruct'.
%
%   In all cases, contact predictions for which there is no data (i.e.,
%   knownmask is false) are ignored in the calculations. Their values in the
%   'precision' ,'tprate', and 'fprate' vectors are set equal to the
%   preceding values if they exist, or to 1 for 'precision' and 'tprate',
%   or 0 for 'fprate'.
%
%   'Check_mapped' is a cell matrix of size equal to the number of
%   alphabets in which each cell has the same format as 'check'. The
%   indices in the fields of check_mapped{i, j} correspond to indices in
%   contacts_pred.mapped{i, j}.
%
%   Options:
%    'threshold' <x>
%       Threshold on the distance between two residues for considering them
%       to be in contact.
%       (default: 6)
%
% See also: GETCONTACTS, GETDISTSTRUCT.

% Tiberiu Tesileanu (2014)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParamValue('threshold', 6, @(x) isscalar(x) && isnumeric(x));

% parse
parser.parse(varargin{:});
params = parser.Results;

nalphas = length(contacts.refseq);
if nalphas ~= length(diststruct.refseq)
    error([mfilename ':alphamissmatch'], 'The predicted contacts and the distance structure should have the same number of alphabets.');
end

% threshold the distance matrix
realcontactmap = (diststruct.distmat < params.threshold);

checked = contacts;

structalphawidths = arrayfun(@(r) length(r.map), diststruct.refseq);
structalpharanges = getalpharanges(structalphawidths);

contactsalphawidths = arrayfun(@(r) length(r.map), contacts.refseq);
contactsalpharanges = getalpharanges(contactsalphawidths);

structidxs = [];
contactsidxs = [];
cidxs = cell(nalphas, 1);
sidxs = cell(nalphas, 1);
for a = 1:nalphas
    [cidxs{a}, sidxs{a}] = findcommonidxs(contacts.refseq(a).map, diststruct.refseq(a).map);
    
    structidxs = [structidxs ; structalpharanges(1, a) + sidxs{a}(:) - 1]; %#ok<AGROW>
    contactsidxs = [contactsidxs ; contactsalpharanges(1, a) + cidxs{a}(:) - 1]; %#ok<AGROW>
end

common_cmap = realcontactmap(structidxs, structidxs);
common_valid_struct = diststruct.validmask(structidxs, structidxs);

% contactsidxs maps indices in common_cmap to raw contacts indices
% now we need the reverse mapping, from raw positions in the contacts
% structure to indices in common_cmap
contactsidxs_inv = zeros(size(contacts.matrix, 1), 1);
contactsidxs_inv(contactsidxs) = 1:length(contactsidxs);

common_raw_pairs = contactsidxs_inv(contacts.rawpairs);

checked.check = check_single_alpha_pair(common_raw_pairs, common_cmap, common_valid_struct, true);

% now work on the mapped contacts
checked.check_mapped = repmat({struct}, nalphas, nalphas);

% in order to reuse check_single_alpha_pair, we need a way to map the
% mapped contacts back to numeric indices... and we need these indices to
% be within the range of indices that are common between the predicted
% contacts and the structure
invmaps = cell(nalphas, 1);
common_alpha_assign = zeros(size(common_cmap, 1), 1);
for a = 1:nalphas
    crtmap = contacts.refseq(a).map;
    invmaps{a} = containers.Map(crtmap, contactsidxs_inv(contactsalpharanges(1, a) + (0:length(crtmap)-1)));
    
    % also need to know which alphabet each index in common_cmap belongs to
    common_alpha_assign(contactsidxs >= contactsalpharanges(1, a) & contactsidxs <= contactsalpharanges(2, a)) = a;
end

for a = 1:nalphas
    alpha_a_mask = (common_alpha_assign == a);
    for b = 1:nalphas
        unmapped_contacts = [...
            cellfun(@(elem) invmaps{a}(elem), contacts.mapped{a, b}(1, :)) ; ...
            cellfun(@(elem) invmaps{b}(elem), contacts.mapped{a, b}(2, :)) ...
        ];
    
        alpha_b_mask = (common_alpha_assign == b);
        
        % now unmapped_contacts have indices in the whole common_cmap
        % matrix... we need to reduce them to the smaller matrix that
        % refers only to alphabets a and b.
        unmapped_contacts_shifted = [...
            arrayfun(@(i) sum(alpha_a_mask(1:i)), unmapped_contacts(1, :)) ; ...
        	arrayfun(@(j) sum(alpha_b_mask(1:j)), unmapped_contacts(2, :)) ...
        ];
        
        checked.check_mapped{a, b} = check_single_alpha_pair(...
            unmapped_contacts_shifted, common_cmap(alpha_a_mask, alpha_b_mask), ...
            common_valid_struct(alpha_a_mask, alpha_b_mask), (a == b));
    end
end

end

function check = check_single_alpha_pair(raw_contacts, cmap, validmap, same)
% CHECK_SINGLE_ALPHA_PAIR Check contacts for a single pair of alphabets.
%   check_struct = CHECK_SINGLE_ALPHA_PAIR(raw_contacts, cmap, validmap, same)
%   checks the contact predictions for a single pair of alphabets,
%   returning a structure with fields 'correctmask', 'knownmask',
%   'precision', 'tprate', and 'fprate' as described for the main function.
%   
%   The evaluation is done using 'cmap', which is the sub-block of the
%   contact map corresponding to the pair of alphabets we're focusing on,
%   and 'validmap', which is the sub-block of the validmask that
%   corresponds to this pair of alphabets. Both of these should be
%   restricted to indices that are common between the contact predictions
%   and the structure. 'Raw_contacts' is a 2 x k matrix of contacts with
%   indices within 'cmap'.
%
%   When 'same' is true, it is assumed that the two alphabets are
%   identical, in which case only elements above the diagonal of cmap and
%   validmap are taken into consideration. Otherwise all the elements are
%   used.

if same
    ncontacts = sum(flatten(triu(cmap & validmap, 1)));
    nnoncontacts = sum(flatten(triu(~cmap & validmap, 1)));
else
    ncontacts = sum(flatten(cmap & validmap));
    nnoncontacts = sum(flatten(~cmap & validmap));
end

% we have no contact predictions for pairs for which there are 0s in
% raw_contacts
check.knownmask = ~any(raw_contacts == 0, 1);
% we also have no contact predictions for pairs for which validmap is false
check.knownmask(check.knownmask) = validmap(...
    sub2ind(size(validmap), raw_contacts(1, check.knownmask), raw_contacts(2, check.knownmask)));

check.correctmask = false(size(check.knownmask));
check.precision = ones(size(check.knownmask));
check.tprate = ones(size(check.knownmask));
check.fprate = zeros(size(check.knownmask));
for i = 1:length(check.correctmask)
    if check.knownmask(i)
        check.correctmask(i) = cmap(raw_contacts(1, i), raw_contacts(2, i));
        
        subcorrect = check.correctmask(1:i);
        subknown = check.knownmask(1:i);
        ntrue = sum(subcorrect);
        
        check.precision(i) = ntrue / sum(subknown);
        if ncontacts > 0
            check.tprate(i) = ntrue / ncontacts;
        end
        
        nfalse = sum(~subcorrect(subknown));
        
        if nnoncontacts > 0
            check.fprate(i) = nfalse / nnoncontacts;
        end
    elseif i > 1
        check.precision(i) = check.precision(i - 1);
        check.tprate(i) = check.tprate(i - 1);
        check.fprate(i) = check.fprate(i - 1);
    end
end

end