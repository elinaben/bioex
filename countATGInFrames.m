
%load the data
T = readtable("yeast_parameters_table_with_diff_5utr.xls");

%use same loop to:
count_ATGS = zeros(1,180);
count_ATGS_rand = zeros(1,180);

number_of_genes = size(T,1)

for i = 1:size(T,1)
    utr5 = char(T{i,"UTR_5"});
    orf = char(T{i,"ORF_1"});

    utr5_rand = char(T{i,"RAND_UTR_5"});
    orf_rand = char(T{i,"RAND_ORF_90"});

    utr5_len_orig = T{i,"UTR5_LEN_ORIG"};
    if isnan(utr5_len_orig)
        continue;
    end
    
    count_ATGS = countATGInString(utr5, orf, count_ATGS);
    count_ATGS_rand = countATGInString(utr5_rand, orf_rand, count_ATGS_rand);


end

figure;
A = [count_ATGS(1:3:end); count_ATGS_rand(1:3:end)]';
bar(log2(A), 'BarWidth', 0.8);
figure;
B = [count_ATGS(2:3:end); count_ATGS_rand(2:3:end)]';
bar(B, 'BarWidth', 0.8);
figure;

C = [count_ATGS(3:3:end); count_ATGS_rand(3:3:end)]';
bar(C, 'BarWidth', 0.5);


for i=1:180
    % define the number of occurrences of "ATG" in the random set and highly expressed genes
    atgCountRandom = count_ATGS(i);
    atgCountRealGenes = count_ATGS_rand(i);
    
    % define the total number of genes in each set
    nRandomSet = 20 * size(T,1);
    nRealSet = size(T,1);

    % calculate the mean and standard deviation of the number of "ATG" occurrences in each set
    meanRandomSet = atgCountRandom / nRandomSet;
    meanRealSet = atgCountRealGenes / nRealSet;
    stdDevRandomSet = sqrt(meanRandomSet * (1 - meanRandomSet) / nRandomSet);
    stdDevHighlyExpressed = sqrt(meanRealSet * (1 - meanRealSet) / nRealSet);
    
    % perform a two-sample t-test to compare the means of the two sets
    [h, pValue] = ttest2([repmat(atgCountRandom, nRandomSet, 1); zeros(nRealSet, 1)], ...
                          [zeros(nRandomSet, 1); repmat(atgCountRealGenes, nRealSet, 1)], ...
                          'Vartype', 'unequal')
end

function count_ATGS = countATGInString(utr5, orf, count_ATGS)
    
    offset_len = 0;
    if length(utr5) < 90
        offset_len = 90 - length(utr5);
    end

    utr5_and_orf = strcat(utr5, orf(1:min(end,90)));
    k = strfind(utr5_and_orf,"ATG");
    
    for j=1:numel(k)
        atg_loc = k(j) + offset_len;
        count_ATGS(atg_loc) = count_ATGS(atg_loc) + 1;
    end
end