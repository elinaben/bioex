
%load the data
T = readtable("yeast_parameters_table.xls");

%remove rows without PA1
colPA1 = rmmissing(T(:,3));
%to get the most abdundant 5% get 95 percentile
P = prctile(table2array(colPA1),95);

%count number of codons in abundant genes - the reference
%convery PA2,PA3 to numeric columns
dict_abundant_codon_freq = dictionary();
for i = 1:size(T,1)
    
    T.PA2num(i) = str2double(char(T.PA2(i)));
    T.PA3num(i) = str2double(char(T.PA3(i)));
    gene_name = char(T{i,"ORF"});
    orf = char(T{i,"ORF_1"});
    
    if T{i,"PA1"} > P
        for j=1:3:size(orf,2)-2
            codon = orf(j:j+2);
            if dict_abundant_codon_freq.isConfigured && dict_abundant_codon_freq.isKey(codon)
                dict_abundant_codon_freq(codon) = dict_abundant_codon_freq(codon) + 1;
            else
                dict_abundant_codon_freq(codon) = 1;
            end
        end
    end
end

%load from file the mapping from codon to amino acid
C = readtable("codons_chart.csv");
dict_codon_to_aa = dictionary();

%for each aa find the most frequent codon 
dict_aa_max_freq_and_sum = containers.Map();
for i = 1:size(C,1)
    codon = char(C{i,1});
    aa = char(C{i,3});

    dict_codon_to_aa(codon) = aa;
    if dict_abundant_codon_freq.isKey(codon)
        C.CODON_FREQ{i} = dict_abundant_codon_freq(codon);
                
        if dict_aa_max_freq_and_sum.isKey(aa) 
            
            aa_max_freq_and_sum = dict_aa_max_freq_and_sum(aa);
            aa_sum = aa_max_freq_and_sum{1} + dict_abundant_codon_freq(codon);
            aa_max = aa_max_freq_and_sum{2};

            if dict_abundant_codon_freq(codon) > aa_max
                aa_max = dict_abundant_codon_freq(codon);
            end

            dict_aa_max_freq_and_sum(aa) = {aa_sum, aa_max};

        else
            aa_max_freq_and_sum = {dict_abundant_codon_freq(codon), dict_abundant_codon_freq(codon)}
            dict_aa_max_freq_and_sum(aa) = aa_max_freq_and_sum;
        end
    else
        C.CODON_FREQ{i} = 0;
    end
    
end

% for each codon its weight according to reference is
% count_in_reference/max_freq_of_codon
dict_codon_weight_and_prob = containers.Map();
codons = keys(dict_codon_to_aa);
for i = 1:size(codons)
    codon = codons(i);
    aa = dict_codon_to_aa(codon);
    if dict_abundant_codon_freq.isKey(codon)
        aa_max_freq_and_sum = dict_aa_max_freq_and_sum(aa);

        weight_and_prob = { dict_abundant_codon_freq(codon) / aa_max_freq_and_sum{2}, dict_abundant_codon_freq(codon) / aa_max_freq_and_sum{1} };
        dict_codon_weight_and_prob(codon) = weight_and_prob ;
    end
end

for i = 1:size(C,1)
    codon = char(C{i,1});
    aa = char(C{i,3});

    weight_prob = dict_codon_weight_and_prob(codon);
    C.CODON_WEIGHT{i} = weight_prob{1};
    C.CODON_FREQ_PERC{i} = weight_prob{2};
end
writetable(C, 'codons_chart_freq.xls')

% for each gene calculate CAI = geomean(codon_weights)
for i = 1:size(T,1)
    gene_name = char(T{i,"ORF"});
    
    orf = char(T{i,"ORF_1"});
    weights = [1 2];
    for j=1:3:size(orf,2)-2
        codon = orf(j:j+2);
        weight = dict_codon_weight_and_prob(codon);
        weights(1+(j-1)/3) = weight{1};
    end

    cai = geomean(weights);
    T.CAI(i) = cai;
end

PA_and_CAI = table2array(rmmissing(T(:,["PA1", "PA2num", "PA3num","CAI"])));
CAI = PA_and_CAI(:,4);
PA_levels = PA_and_CAI(:,1:3);
reg_corr = corr(PA_levels,CAI)
spearman_corr = corr(PA_levels, CAI, "type", "Spearman")
scatter(log(PA_levels(:,1)), CAI)

