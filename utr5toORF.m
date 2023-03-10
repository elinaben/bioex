filetextOrfCodingAll = fileread('orf_coding_all.fasta');
T = readtable("yeast_parameters_table.xls");
T_UTR5 = readtable("TableS4.xls");

%load from file the mapping from codon to amino acid with probability of
%frequent codons
[aa_to_codon_freq, codon_aa] = get_codon_probabilities();

chr_file = 0;
nucleotides = '';
missing_count = 0;
for i = 1:2%size(T,1)
    
    gene_just_name = char(T{i,"ORF"});
    gene_name = ">" + gene_just_name + " ";
    gene_location_start = strfind(filetextOrfCodingAll,gene_name);
    header = filetextOrfCodingAll(gene_location_start:gene_location_start+100);
    T.LOC{i} = header;
    if (T.LOC{i} == "")
        missing_count = missing_count + 1
        disp(gene_name)
        continue
    end 
    [nucleotides, chr_file] = read_nucleaotides(header, nucleotides, chr_file);
    start_end_array = getStartAndEndLocations(T.LOC{i});

    orf_2 = '';
    %concatenate all the ORF's
    for j = 1:size(start_end_array,1)
        start_i = start_end_array(j,1);
        end_i = start_end_array(j,2);   

        if start_i < end_i
            orf_2 = [orf_2, nucleotides(start_i:end_i)];   
        else
            sense = getSenseStrand(nucleotides(end_i:start_i));
            orf_2 = [sense, orf_2];
        end
    end

    % find utr5' length and concatenate to ORF
    row_idx = strcmp(T_UTR5.Name, gene_just_name);

    % if we found the gene in the file
    if sum(row_idx) > 0
        x5__UTR_length = T_UTR5{row_idx,"x5__UTR_length"};
        %if the length doesn't appear use 90 as default
        if isnan(x5__UTR_length)
            x5__UTR_length = 90;
        end

        %don't take more than 90nt
        if x5__UTR_length > 90
            x5__UTR_length = 90;
        end
        
        %if the gene is on this strand take the utr before the orf
        %otherwise we read the part after and then get the sense from the
        %anti sense
        if start_end_array(1,1) < start_end_array(1,2)
            if start_i - x5__UTR_length - 1 > 0
                utr_5 = nucleotides(start_i-x5__UTR_length:start_i-1);
            else 
                disp("missing full utr: " + gene_name)
            end
        else
            start_i = start_end_array(end,1);
            utr_5 = getSenseStrand(nucleotides(start_i+1:start_i+x5__UTR_length));
        end
    end
    T.UTR_5{i} = utr_5;
    T.UTR5_LEN{i} = x5__UTR_length;
    T.UTR5_LEN_ORIG{i} = T_UTR5{row_idx,"x5__UTR_length"};

    T.RAND_UTR_5{i} = get_random_utr(utr_5);
    T.RAND_ORF{i} = get_random_orf(orf_2(1:90), aa_to_codon_freq, codon_aa);
    T.RAND_ORF_2{i} = get_random_orf(orf_2(1:90), get_orf_probabilities(orf_2, codon_aa), codon_aa);

end

% calculate random 20 ORFs, count ATGs and create an avaverage 
% how do I calculate the p-value?


%writetable(T, 'yeast_parameters_table_with_diff_5utr.xls')

function [aa_to_codon_freq, codon_aa] = get_codon_probabilities()
    %load from file the mapping from codon to amino acid
    C = readtable("codons_chart_freq.xls");
    
    aa_to_codon_freq = containers.Map();
    codon_aa = containers.Map();
    for i = 1:size(C,1)
        codon = char(C{i,1});
        aa = char(C{i,3});
    
        codon_aa(codon) = aa;
        if ~aa_to_codon_freq.isKey(aa)
            index = cellfun(@(x) isequal(x, aa), C.aa);
            arrCodons = C.codon(index);
            arrProb = C.CODON_FREQ_PERC(index);
    
            aa_to_codon_freq(aa) = {arrCodons, arrProb};
        end
    end

end

function aa_to_codon_freq = get_orf_probabilities(orf, codon_aa)

    
    aa_to_codon_freq = containers.Map();

    codons_count = containers.Map();
    aa_count = containers.Map();

    for j=1:3:size(orf,2)-2
        codon = orf(j:j+2);
        aa = codon_aa(codon);
        if ~codons_count.isKey(codon)
            codons_count(codon) = 1;
        else
            codons_count(codon) = codons_count(codon) + 1;
        end

        if ~aa_count.isKey(aa)
            aa_count(aa) = 1;
        else
            aa_count(aa) = aa_count(aa) + 1;
        end
    end

    codons = keys(codons_count);
    for i=1:numel(codons)
        codon = char(codons(i));
        aa = codon_aa(codon);
        if ~aa_to_codon_freq.isKey(aa)
            arrCodons = [codon];
            arrProb = [codons_count(codon)/aa_count(aa)];
    
            aa_to_codon_freq(aa) = {arrCodons, arrProb};
        else
            arrCodons = [arrCodons,""+codon]
            arrProb = [arrProb, codons_count(codon)/aa_count(aa)]
    
            aa_to_codon_freq(aa) = {arrCodons, arrProb};
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

function [nucleotides_out, curr_chr_file] = read_nucleaotides(header, nucleotides, prev_file_number)

    chr_start = strfind(header, "Chr");
    chr_end = strfind(header, "from");
    chr_num = header(chr_start+4:chr_end-2);
    switch chr_num
        case 'I'
            curr_chr_file = 1;
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr01.fsa';
        case 'II'
            curr_chr_file = 2;
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr02.fsa';
        case 'III'
            curr_chr_file = 3;
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr03.fsa';
        case 'IV'
            curr_chr_file = 4;     
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr04.fsa';
        case 'V'
            curr_chr_file = 5;     
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr05.fsa';
        case 'VI'
            curr_chr_file = 6;     
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr06.fsa';
        case 'VII'
            curr_chr_file = 7;     
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr07.fsa';
        case 'VIII'
            curr_chr_file = 8;
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr08.fsa';
        case 'IX'
            curr_chr_file = 9;    
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr09.fsa';
        case 'X'
            curr_chr_file = 10;    
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr10.fsa';
        case 'XI'
            curr_chr_file = 11;    
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr11.fsa';
        case 'XII'
            curr_chr_file = 12;    
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr12.fsa';
        case 'XIII'
            curr_chr_file = 13; 
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr13.fsa';
        case 'XIV'
            curr_chr_file = 14;    
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr14.fsa';
        case 'XV'
            curr_chr_file = 15;    
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr15.fsa';
        case 'XVI'
            curr_chr_file = 16;    
            curr_chr_file_name = 'http://sgd-archive.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr16.fsa';
           
    end
        
    if prev_file_number ~= curr_chr_file
        
        %open the chromosome file
        options = weboptions("ContentType",'text');
        file_content = webread(curr_chr_file_name, options);
                    
        nucleotides_out = '';
        lines = strsplit(file_content, '\n');
        for j = 1:length(lines)
            if isempty(lines{j}) || lines{j}(1) == '>'
                continue;
            end
            nucleotides_out = [nucleotides_out, lines{j}];
        end
    else
        nucleotides_out = nucleotides;
    end
    
end

function start_end_array = getStartAndEndLocations(header)

    gene_location_start = strfind(header,"from");
    gene_location_end = strfind(header,"Genome");
    header = header(gene_location_start + 4:gene_location_end - 3);

    range_cell = split(header, ",");
    n = length(range_cell);
    start_end_array = zeros(n,2);
    for i = 1:n
        range = range_cell{i};
        numbers = split(range, '-');
        start_end_array(i,:) = str2double(numbers);
    end

end


function complementary_strand = getSenseStrand(anti_sense_strand)
    complementary_strand = "";
    for i = 1:length(anti_sense_strand)
        switch anti_sense_strand(i)
            case 'A'
                complementary_strand(i) = 'T';
            case 'T'
                complementary_strand(i) = 'A';
            case 'C'
                complementary_strand(i) = 'G';
            case 'G'
                complementary_strand(i) = 'C';
        end
    end
    complementary_strand = strjoin(cellstr(complementary_strand),'');
    complementary_strand = reverse(complementary_strand);
end