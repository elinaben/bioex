T = readtable("yeast_parameters_table_with_diff_5utr.xls");

[aa_to_codon, codon_to_aa] = get_aa_to_codon_map()

randUTRandORFTable = table('Size',[0 5],'VariableTypes',{})

for rand_count = 1:20
    for i = 1:size(T,1)
        
        orf = char(T{i,"ORF_1"});
        utr5 = char(T{i,"UTR_5"});
        
        utr5_rand = char(T{i,"RAND_UTR_5"});
        orf_rand = char(T{i,"RAND_ORF_90"});
    
        utr5_len_orig = T{i,"UTR5_LEN_ORIG"};
        if isnan(utr5_len_orig)
            continue;
        end
        
        aa_keys = keys(aa_to_codon);
        codon_count = containers.Map();
        rand_orf = orf;
        for j = 1:numel(aa_keys)
            codons = aa_to_codon(char(aa_keys(j)));
            idxes = [];
            for k=1:numel(codons)
                codon_idxes = strfind(orf, codons(k));
                codon_count(char(codons(k))) = size(codon_idxes);
                idxes = [idxes, codon_idxes];
            end
    
            %generate a random permutation of the indices
            perm_idxes = randperm(length(idxes));
            idxes_after_perm = idxes(perm_idxes);
    
            for k=1:numel(idxes_after_perm)
                pre_idx = idxes(k);
                new_idx = idxes_after_perm(k);
                if strcmp(orf(pre_idx:pre_idx+2), orf(new_idx:new_idx+2)) == 0
                    rand_orf(new_idx:new_idx+2) = orf(pre_idx:pre_idx+2);
                end
            end
    
        end
    
        T.RAND_ORF_90{i} = rand_orf;
    end




end

% calculate random 20 ORFs, count ATGs and create an avaverage 
% how do I calculate the p-value?

%writetable(T, 'yeast_parameters_table_with_diff_5utr.xls')

function [aa_to_codon, codon_to_aa] = get_aa_to_codon_map()
    %load from file the mapping from codon to amino acid
    C = readtable("codons_chart_freq.xls");
    
    aa_to_codon = containers.Map();
    codon_to_aa = containers.Map();
    for i = 1:size(C,1)
        codon = char(C{i,1});
        aa = char(C{i,3});
    
        codon_to_aa(codon) = aa;
        if ~aa_to_codon.isKey(aa)
            index = cellfun(@(x) isequal(x, aa), C.aa);
            arrCodons = C.codon(index);
            aa_to_codon(aa) = arrCodons;
        end
    end
    
end

function rand_orf = get_random_orf(orf, aa_to_codon_freq, codon_aa)
    
    % for each amino acid create a vector that contains all the locations
    % of that a.a. and then create a random sample with the size according
    % to number of this a.a. in original orf and replace them with
    % the random codons for that a.a.
    aa_to_idx = containers.Map();
    for j=1:3:size(orf,2)-2
        codon = orf(j:j+2);
        aa = codon_aa(codon);
        if ~aa_to_idx.isKey(aa)
            aa_to_idx(aa) = [];
        end

        aa_to_idx(aa) = [aa_to_idx(aa), j];
    end
    rand_orf = orf;
    aa_keys = keys(aa_to_codon_freq);
    for j=1:length(aa_keys)

        aa2 = char(aa_keys(j));
        if ~aa_to_idx.isKey(aa2)
            continue
        end
        aa_indices = aa_to_idx(aa2);
        codons_count = length(aa_indices);
        codons_and_freq = aa_to_codon_freq(aa2);

        if codons_and_freq{2} == 1
            continue
        end

        sample = randsample(codons_and_freq{1}, codons_count, true, codons_and_freq{2});

        for k=1:codons_count
            idx = aa_indices(k);
            rand_orf = char(strcat(rand_orf(1:idx-1), sample(k), rand_orf(idx+3:end)));
        end
    end
end


function ran_utr = get_random_utr(utr)
    n = length(utr);
    %generate a random permutation of the indices
    idx = randperm(n);
    ran_utr = utr(idx);
end