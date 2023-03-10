
%load the data
T = readtable("yeast_parameters_table_with_5utr.xls");

%remove rows without PA1
colPA1 = rmmissing(T(:,3));
%to get the most abdundant 5% get 95 percentile
P = prctile(table2array(colPA1),95);

%use same loop to:
% 1 - count per gene how many times we saw each codon
% 2 - count number of codons in abundant genes - the reference
dict_abundant_codon_freq = dictionary();
dict_gene_codon_freq = dictionary();
dict_gene_to_num_of_codons = dictionary();
for i = 1:size(T,1)
    gene_name = char(T{i,"ORF"});
    orf = char(T{i,"ORF_1"});
    [gene_codon_freq, gene_number_of_codons] = CountCodonsInGene(orf);
    dict_gene_codon_freq(gene_name) = gene_codon_freq;
    dict_gene_to_num_of_codons(gene_name) = gene_number_of_codons;
    %use only the top 5 percent as a reference to calculate freq
    %if the gene is in reference add all the counters of it's codons to the
    %counter of reference codons
    if T{i,"PA1"} > P
        codons = keys(gene_codon_freq);
        for codon_i = 1:size(codons)
            codon = codons(codon_i);
            if dict_abundant_codon_freq.isConfigured && dict_abundant_codon_freq.isKey(codon)
                dict_abundant_codon_freq(codon) = dict_abundant_codon_freq(codon) + gene_codon_freq(codon);
            else
                dict_abundant_codon_freq(codon) = gene_codon_freq(codon);
            end
        end
    end
end

%load from file the mapping from codon to amino acid
C = readtable("codons_chart.csv");

dict_codon_to_aa = dictionary();

%for each aa find the most frequent codon and store it's count
%in a map from aa to the count of most frequent codon
dict_aa_max_freq = dictionary();
for i = 1:size(C,1)
    codon = char(C{i,1});
    aa = char(C{i,3});

    dict_codon_to_aa(codon) = aa;
    if dict_abundant_codon_freq.isKey(codon)
        if dict_aa_max_freq.isConfigured && dict_aa_max_freq.isKey(aa) 
            if dict_abundant_codon_freq(codon) > dict_aa_max_freq(aa)
                dict_aa_max_freq(aa) = dict_abundant_codon_freq(codon);
            end
        else
            dict_aa_max_freq(aa) = dict_abundant_codon_freq(codon);
        end
    end
end

% for each codon it's weight according to reference is
% count_in_reference/max_freq_of_codon
dict_codon_weight = dictionary();
codons = keys(dict_codon_to_aa);
for i = 1:size(codons)
    codon = codons(i);
    aa = dict_codon_to_aa(codon);
    if dict_abundant_codon_freq.isKey(codon)
        dict_codon_weight(codon) = dict_abundant_codon_freq(codon) / dict_aa_max_freq(aa);
    end
end


for i = 1:size(T,1)
    gene_name = char(T{i,"ORF"});
    orf = char(T{i,"ORF_1"});
    
    for j=1:3:size(orf,2)-2
        codon = orf(j:j+2);
        
    end
end

% for each gene calculate CAI = geomean(codon_weights)
% instead of multipling weights - which can be inaccurate because the
% numbers become very small use the formula:
% single_codon_weight^(num_of_this_codon/number_of_total_codons)
gene_names = keys(dict_gene_codon_freq);
dict_gene_to_cai = dictionary();
for gene_name_i = 1:size(gene_names)
    gene_name = gene_names(gene_name_i)
    gene_codon_freq = dict_gene_codon_freq(gene_name)
    
    total_cai = 1;
    codons = keys(gene_codon_freq);
    for i = 1:size(codons)
       codon = codons(i);
       power_codon = gene_codon_freq(codon)/dict_gene_to_num_of_codons(gene_name)
       weight = dict_codon_weight(codon)
       cai = power(weight,power_codon)
       total_cai = total_cai*cai
    end
    dict_gene_to_cai(gene_name,cai)
end

function [dict_codons_count,codons_count] = CountCodonsInGene(orf)
    dict_codons_count = dictionary();
    codons_count = 0;
    for j=1:3:size(orf,2)-2
        codon = orf(j:j+2);
        codons_count = codons_count + 1;
        if dict_codons_count.isConfigured && dict_codons_count.isKey(codon)
            dict_codons_count(codon) = dict_codons_count(codon) + 1;
        else
            dict_codons_count(codon) = 1;
        end
    end
end