

T = readtable("yeast_parameters_table_with_diff_5utr.xls");
[aa_to_codon, codon_to_aa] = get_aa_to_codon_map();

SLIDING_WINDOW_SIZE = 40;
max_length = 80;
ORF_SIZE = 550;

folding_energy_array = zeros(size(T,1), SLIDING_WINDOW_SIZE + ORF_SIZE);
for i = 1:size(T,1)
    disp(i);
    orf = char(T{i,"ORF_1"});
	orf_rand_sample = getRandORF(orf, aa_to_codon);
    utr5 = char(T{i,"UTR_5"});
    utr5_len_orig = T{i,"UTR5_LEN_ORIG"};
    offset = 1;
    if length(utr5) < max_length
        offset = max_length - length(utr5) + 1;
    elseif length(utr5) > max_length
        %concatenate utr5 to take only the last max_length chars
        utr5 = utr5(length(utr5)-max_length+1:end);
    end
    utr5_rand_sample = getRandUTR(utr5);

    full_sequence = strcat(utr5_rand_sample,orf_rand_sample(1:min(end,ORF_SIZE)));
    disp(full_sequence);
    folding_energy_result = runMFold(full_sequence, SLIDING_WINDOW_SIZE);
    folding_energy_array(i,offset:(offset + length(folding_energy_result) - 1)) = folding_energy_result;
    
    if mod(i,10) == 0
        disp("writing...");
        csvwrite('folding_energy_random.csv',folding_energy_array(1:i,:));
    end
end

mean_values = mean(folding_energy_array);
std_values = std(folding_energy_array);

figure;
plot(mean_values, "LineWidth", 2);
figure;
plot(std_values);

csvwrite("folding_energy.csv",folding_energy_array);

function num = extract_last_parenthesis_number(str)
    % Extract all numbers within parentheses
    matches = regexp(str, '\((\s*-?\d+(\.\d+)?)\)', 'tokens');
    num = 0.0;
    % Check if any matches are found
    if isempty(matches)
        disp('No number found in the last parenthesis: ' + str);
        return;
    end
    
    % Extract the last match (last number in parentheses)
    last_match = matches{end};
    num = str2double(last_match{1});
end

function folding_energy_array = runMFold(orf, sliding_window_size)

    folding_energy_array = zeros(1, length(orf));
    
    for j=1:length(orf)
        sequence = orf(j:min(j+sliding_window_size,length(orf)));
        command = sprintf('echo "%s" | RNAfold ', sequence);
        [result, folding_energy] = system(command);
        
        
        num = extract_last_parenthesis_number(folding_energy);
        folding_energy_array(j) = num;
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

function folding_energy_array = runMFold_redirect(orf, sliding_window_size)

    folding_energy_array = zeros(1, length(orf));
    
    for j=1:length(orf)
        sequence = orf(j:min(j+sliding_window_size,length(orf)));
        command = sprintf('echo "%s" | RNAfold > myFile_%d', sequence, j);
        [result, folding_energy] = system(command);
        
        
        %num = extract_last_parenthesis_number(folding_energy);
        folding_energy_array(j) = j;
    end
end