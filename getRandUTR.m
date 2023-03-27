function ran_utr = getRandUTR(utr)
    n = length(utr);
    %generate a random permutation of the indices
    idx = randperm(n);
    ran_utr = utr(idx);
end


