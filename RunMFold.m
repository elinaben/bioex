
T = readtable("yeast_parameters_table_with_diff_5utr.xls");
SLIDING_WINDOW_SIZE = 40;
MAX_UTR5_LEN = 23;
ORF_SIZE = 550;

folding_energy_array = zeros(size(T,1), MAX_UTR5_LEN + ORF_SIZE);
%T_results = readmatrix("folding_energy_30_03.csv");
%for i = 1:size(T_results,1)%for i = size(T,1)-1:-1:630
%    folding_energy_array(i,:) = T_results(i,:);
%end
%size(T_results,1)+1:size(T,1)
for i = 1:size(T,1)
    disp(i);
    orf = char(T{i,"ORF_1"});
    utr5 = char(T{i,"UTR_5"});
    utr5_len_orig = T{i,"UTR5_LEN_ORIG"};
    offset = 1;
    if length(utr5) < MAX_UTR5_LEN
        offset = MAX_UTR5_LEN - length(utr5) + 1;
    elseif length(utr5) > MAX_UTR5_LEN
        %concatenate utr5 to take only the last max_length chars
        utr5 = utr5(length(utr5)-MAX_UTR5_LEN+1:end);
    end
    
    full_sequence = strcat(utr5,orf(1:min(end,ORF_SIZE)));
    disp(full_sequence);
    folding_energy_result = runMFold(full_sequence, SLIDING_WINDOW_SIZE);
    folding_energy_array(i,offset:(offset + length(folding_energy_result) - 1)) = folding_energy_result;
    
    if mod(i,10) == 0
        disp("writing...");
        file_name = sprintf('folding_energy_bioinformatics.csv', i);
        csvwrite(file_name,folding_energy_array(1:i,:));%size(T,1)-1,:));
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
        %command = sprintf('echo "%s" | RNAfold ', sequence);
        %[result, folding_energy] = system(command);
        [RNAbracket, Energy] = rnafold(sequence);
                
        folding_energy_array(j) = Energy;
    end
end
