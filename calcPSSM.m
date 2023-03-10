%load the data
T = readtable("yeast_parameters_table_with_5utr.xls");

%remove rows without PA1
colPA1 = rmmissing(T(:,3));
%to get the most abdundant 4% get 96 percentile
P = prctile(table2array(colPA1),96)
T_abundant = T(T.PA1 >= P, :);

% Initialize the numeric matrix with NaN values
numseqs = nan(height(T_abundant), 4);
cnt_nt_pos = zeros(4,12);

for i = 1:size(T_abundant,1)
    % Extract the 5' UTR sequence
    orf_utr5 = char(T_abundant{i,"UTR_5"});
    orf_utr5 = orf_utr5(86:98);

    for j = 1:12
        if orf_utr5(j) == 'C'
            cnt_nt_pos(1,j) = cnt_nt_pos(1,j) + 1;
        elseif orf_utr5(j) == 'A'
            cnt_nt_pos(2,j) = cnt_nt_pos(2,j) + 1;
        elseif orf_utr5(j) == 'T'
            cnt_nt_pos(3,j) = cnt_nt_pos(3,j) + 1;
        elseif orf_utr5(j) == 'G'
            cnt_nt_pos(4,j) = cnt_nt_pos(4,j) + 1;
            
        end
    end
    
end



cnt_nt_pos = cnt_nt_pos';
row_sums = sum(cnt_nt_pos,2);
pssm = cnt_nt_pos ./ row_sums;

colorVector = ["blue", "yellow", "magenta", "green"];
b = bar(pssm, 'stacked', 'BarWidth', 1);
for i = 1:numel(b)
    b(i).FaceColor = colorVector(i);
end
figure; 

entropy = -sum(pssm .* log2(pssm), 2);

cnt_nt_pos_2 = zeros(4,12);
for i = 1:size(T, 1)
    
    T.PA1num(i) = str2double(char(T.PA1(i)));
    T.PA2num(i) = str2double(char(T.PA2(i)));
    T.PA3num(i) = str2double(char(T.PA3(i)));
    
    % Extract the 5' UTR sequence
    orf_utr5 = char(T{i,"UTR_5"});
    if length(orf_utr5) < 57
        continue;
    end
    orf_utr5 = orf_utr5(86:98);

    % calculate the context score by using the probability for this nt from
    % pssm
    context_score = 0;
    for j = 1:12
        col = 1; 
        if orf_utr5(j) == 'C'
            col = 1;
            cnt_nt_pos_2(1,j) = cnt_nt_pos_2(1,j) + 1;
        elseif orf_utr5(j) == 'A'
            col = 2;
            cnt_nt_pos_2(2,j) = cnt_nt_pos_2(2,j) + 1;
        elseif orf_utr5(j) == 'T'
            col = 3;
            cnt_nt_pos_2(3,j) = cnt_nt_pos_2(3,j) + 1;
        elseif orf_utr5(j) == 'G'
            col = 4;
            cnt_nt_pos_2(4,j) = cnt_nt_pos_2(4,j) + 1;
        end

        log_val = log2(pssm(j,col));
        context_score = context_score + log_val;
    end
    
    T.CONTEXT_SCORE(i) = context_score;

    T.CONTEXT_SCORE_EXP(i) = exp(context_score);
end


cnt_nt_pos_2 = cnt_nt_pos_2';
row_sums_2 = sum(cnt_nt_pos_2,2);
pssm_2 = cnt_nt_pos_2 ./ row_sums_2;
b2 = bar(pssm_2, 'stacked', 'BarWidth', 1);
for i = 1:numel(b2)
    b2(i).FaceColor = colorVector(i);
end
figure; 


PA_and_CNXT = table2array(rmmissing(T(:,["PA1", "PA2num", "PA3num","CONTEXT_SCORE"])));
CONTEXT_SCORE = PA_and_CNXT(:,4);
PA_levels = PA_and_CNXT(:,1:3);
reg_corr = corr(PA_levels,CONTEXT_SCORE)
spearman_corr = corr(PA_levels, CONTEXT_SCORE, "type", "Spearman")
%scatter(log(PA_levels(:,1)), CONTEXT_SCORE)
