function rand_orf = getRandomORF(orf, aa_to_codon)
    
    aa_keys = keys(aa_to_codon);
    rand_orf = orf;
    for j = 1:numel(aa_keys)
        codons = aa_to_codon(char(aa_keys(j)));
        idxes = [];
        for k=1:numel(codons)
            codon_idxes = strfind(orf, codons(k));
            codon_idxes = codon_idxes(mod(codon_idxes,3) == 1);
            idxes = [idxes, codon_idxes];
        end

        %generate a random permutation of the indices
        perm_idxes = randperm(length(idxes));
        idxes_after_perm = idxes(perm_idxes);

        for k=1:numel(idxes_after_perm)
            pre_idx = idxes(k);
            new_idx = idxes_after_perm(k);
            pre_codon = orf(pre_idx:pre_idx+2);
            new_codon = orf(new_idx:new_idx+2);
            if strcmp(pre_codon, new_codon) == 0
                rand_orf(new_idx:new_idx+2) = new_codon;
            end
        end
    end
end



