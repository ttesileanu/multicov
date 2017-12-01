function [idxs1, idxs2] = findcommonidxs(ref1, ref2)
% FINDCOMMONIDXS Find places in two arrays where they contain elements
% that are common to both.
%   [idxs1, idxs2] = FINDCOMMONIDXS(ref1, ref2) returns two arrays of
%   indices idxs1 and idxs2 such that ref1(idxs1(i)) = ref2(idxs2(i)) for
%   all i. The arrays ref1 and ref2 may be numeric arrays, or cell arrays
%   of strings.

% Tiberiu Tesileanu (2013-2014)

% no match between strings and numbers
if (isnumeric(ref1) && iscell(ref2)) || (isnumeric(ref2) && iscell(ref1))
    idxs1 = [];
    idxs2 = [];
else
    common = intersect(ref1, ref2);
    if isnumeric(common)
        % 0 does not actually refer to any position in the reference sequence
        common = setdiff(common, 0);
    end
  
    [mask1, loc1] = ismember(ref1, common);
    [mask2, loc2] = ismember(ref2, common);
    
    idxs1 = zeros(size(common));
    idxs2 = zeros(size(common));
    
    idxs1(loc1(loc1 ~= 0)) = find(mask1);
    idxs2(loc2(loc2 ~= 0)) = find(mask2);
end

end