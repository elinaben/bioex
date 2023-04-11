T = readtable("yeast_parameters_table.xls");
[aa_to_codon, codon_to_aa] = get_aa_to_codon_map();

% Extract the column you want to check
%col = "PA1";

% Find the rows with non-double values
%non_double_rows = isnan(T.(col));

% Remove the non-double rows from the table
%T(non_double_rows,:) = [];

codon_to_tai = get_codon_to_tAI();
CODON_LEN = 200;

array_gene_tAI = zeros(size(T,1), CODON_LEN);
array_gene_tAI_random = zeros(size(T,1), CODON_LEN);
            
for i = 1:size(T,1)
        
    disp(i);
    orf = char(T{i,"ORF_1"});
    for j=4:3:min(size(orf,2)-2,CODON_LEN * 3)
        codon = orf(j:j+2);
        array_gene_tAI(i,(j+2)/3 - 1) = codon_to_tai(codon);
    end
    
    array_gene_tAI_random_temp = zeros(20, CODON_LEN);
    for k=1:1:100    
        orf_rand_sample = getRandORF(orf, aa_to_codon);
        for j=4:3:min(size(orf,2)-2,CODON_LEN * 3)
            codon = orf_rand_sample(j:j+2);
            array_gene_tAI_random_temp(k,(j+2)/3 - 1) = codon_to_tai(codon);
        end
    end

    array_gene_tAI_random_temp(array_gene_tAI_random_temp == 0) = NaN;
    array_gene_tAI_random_mean_vals = mean(array_gene_tAI_random_temp, 'omitnan');  
    array_gene_tAI_random(i,:) = array_gene_tAI_random_mean_vals;
end


csvwrite('array_gene_tAI.csv',array_gene_tAI);
csvwrite('array_gene_tAI_random.csv',array_gene_tAI_random);
        
array_gene_tAI(array_gene_tAI == 0) = NaN;
mean_vals = mean(array_gene_tAI, 'omitnan');

array_gene_tAI_random(array_gene_tAI_random == 0) = NaN;
mean_vals_rand = mean(array_gene_tAI_random, 'omitnan');

figure;
%x = linspace(UTR_SIZE,length(mean_vals),length(mean_vals));
plot(mean_vals(1:200), "LineWidth", 1.5, "Color", "black");
hold on;
plot(mean_vals_rand(1:200), "LineWidth", 1.5, "Color", "blue");

function [codon_to_tai] = get_codon_to_tAI()
    %load from file the mapping from codon to amino acid
    C = readtable("yeast_codon_tAI.csv");
    
    codon_to_tai = containers.Map();
    for i = 1:size(C,1)
        codon = char(C{i,1});
        tai = double(C{i,2});
        codon_to_tai(codon) = tai;
    end
    
end

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