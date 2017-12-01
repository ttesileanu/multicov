function test_loadfasta
% TEST_LOADFASTA Test loadfasta.m.

% test loading for several alphabets
lf_singletest('protein', '-ACDEFGHIKLMNPQRSTVWY');
lf_singletest('dna', '-ACGT');
lf_singletest('rna', '-ACGU');

% test keepmask
lf_singletest('protein', '-ACDEFGHIKLMNPQRSTVWY', 'keepmask_v');
lf_singletest('dna', '-ACGT', 'keepmask_f');
lf_singletest('rna', '-ACGU', 'matclean');
lf_singletest('protein', '-ACDEFGHIKLMNPQRSTVWY', 'postprocess');

% test gaps
lf_singletest('gapprotein', '-ACDEFGHIKLMNPQRSTVWY');

end