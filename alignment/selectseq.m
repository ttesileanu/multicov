function subalign = selectseq(alignment, idxs)
% SELECTSEQ Select a subset of the sequences in an alignment.
%   subalign = SELECTSEQ(alignment, idxs) returns a subalignment of
%   'alignment' containing the sequences indexed by 'idxs'. Note that
%   'idxs' can also be a logical array identifying the sequences to keep.
%   This can be either a character or a binary alignment.

% Tiberiu Tesileanu (2014)

subalign = alignment;
subalign.data = alignment.data(idxs, :);
subalign.seqw = alignment.seqw(idxs);
subalign.annotations = alignment.annotations(idxs);

end