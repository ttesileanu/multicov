function test_fishertest
% TEST_FISHERTEST Test fishertest.m.

% 1 tail
timer = tic;
ntests = 80;
nmax = 32;
correct = true;
for i = 1:ntests
    matrix = randi([0 nmax], [2 2]);
    p_rel_p = fishertest(matrix, 'tails', 1);
    p_rel_m = fishertest(matrix, 'tails', -1);
    
    p_abs_p = fishertest(matrix, 'tails', 1, 'type', 'abs');
    p_abs_m = fishertest(matrix, 'tails', -1, 'type', 'abs');
    
    matrix0 = matrix;
    
    % recalculate p-values
    step = [1 -1 ; -1 1];
    nsteps = min(matrix(step < 0));
    my_p_abs_p = 0;
    v1 = [1:sum(matrix(1, :)) 1:sum(matrix(2, :)) 1:sum(matrix(:, 1)) 1:sum(matrix(:, 2))];
    for j = 1:nsteps+1
        v2 = [1:matrix(1, 1) 1:matrix(1, 2) 1:matrix(2, 1) 1:matrix(2, 2) 1:sum(matrix(:))];
        crt_prob = prod(v1./v2);
        my_p_abs_p = my_p_abs_p + crt_prob;
        matrix = matrix + step;
    end
    
    step = -step;
    matrix = matrix0;
    nsteps = min(matrix(step < 0));
    my_p_abs_m = 0;
    v1 = [1:sum(matrix(1, :)) 1:sum(matrix(2, :)) 1:sum(matrix(:, 1)) 1:sum(matrix(:, 2))];
    for j = 1:nsteps+1
        v2 = [1:matrix(1, 1) 1:matrix(1, 2) 1:matrix(2, 1) 1:matrix(2, 2) 1:sum(matrix(:))];
        crt_prob = prod(v1./v2);
        my_p_abs_m = my_p_abs_m + crt_prob;
        matrix = matrix + step;
    end
    
    if matrix0(1, 1)*matrix0(2, 2) < matrix0(1, 2)*matrix0(2, 1)
        my_p_rel_p = my_p_abs_m;
        my_p_rel_m = my_p_abs_p;
    else
        my_p_rel_p = my_p_abs_p;
        my_p_rel_m = my_p_abs_m;
    end
    
    if norm([p_rel_p p_rel_m p_abs_p p_abs_m]./[my_p_rel_p my_p_rel_m my_p_abs_p my_p_abs_m] - 1) > 1e-12
        correct = false;
        break;
    end
end
utexpect(correct, 'fishertest 1-tailed', timer);

% 2 tails
timer = tic;
ntests = 80;
nmax = 32;
correct = true;
for i = 1:ntests
    matrix = randi([0 nmax], [2 2]);
    p = fishertest(matrix, 'tails', 2);
    
    % recalculate the p-value
    v1 = [1:sum(matrix(1, :)) 1:sum(matrix(2, :)) 1:sum(matrix(:, 1)) 1:sum(matrix(:, 2))];
    v2 = [1:matrix(1, 1) 1:matrix(1, 2) 1:matrix(2, 1) 1:matrix(2, 2) 1:sum(matrix(:))];
    crt_p0 = prod(v1./v2);
    step = [1 -1; -1 1];
    nsteps = min(matrix(step < 0));
    matrix = matrix + nsteps*step;
    step = -step;
    nsteps = min(matrix(step < 0));
    my_p = 0;
    for j = 0:nsteps
        v2 = [1:matrix(1, 1) 1:matrix(1, 2) 1:matrix(2, 1) 1:matrix(2, 2) 1:sum(matrix(:))];
        crt_prob = prod(v1./v2);
        if crt_prob <= crt_p0 + eps
            my_p = my_p + crt_prob;
        end
        matrix = matrix + step;
    end
    
    if abs(p/my_p - 1) > 1e-12
        correct = false;
        break;
    end
end
utexpect(correct, 'fishertest 2-tailed', timer);

end