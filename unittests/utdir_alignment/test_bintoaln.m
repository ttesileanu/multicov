function test_bintoaln
% TEST_BINTOALN Test bintoaln.m.

% reproducible arbitrariness
old_state = rng;

rng(13);

alignment1 = alngenrandom(32, 10, 'rna');
binalign1 = alntobin(alignment1);
recoveredalign1 = bintoaln(binalign1);
utexpect(alncheck(recoveredalign1) && isequal(recoveredalign1.data, alignment1.data), ...
    'bintoaln rna');

alignment2 = alngenrandom(32, 12, 'protein');
alignment = alnadd(alignment1, alignment2);
binalign2 = alntobin(alignment2);
binalign = alntobin(alignment);
recoveredalign2 = bintoaln(binalign2);
recoveredalign = bintoaln(binalign);
utexpect(alncheck(recoveredalign2) && alncheck(recoveredalign) && ...
    isequal(recoveredalign.data, [recoveredalign1.data recoveredalign2.data]) && ...
    isequal(recoveredalign2.data, alignment2.data), 'bintoaln multi-alphabet');

rng(old_state);

end
