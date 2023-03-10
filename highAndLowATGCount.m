T = readtable("yeast_parameters_table_with_diff_5utr.xls");

%remove rows without PA1
colPA1 = rmmissing(T(:,3));
%to get the most abdundant 4% get 96 percentile
P_high = prctile(table2array(colPA1),85)
P_low = prctile(table2array(colPA1),15)

atg_counter_len = 180;
utr_orf_len = 90;

count_ATGS_high_PA = zeros(1,atg_counter_len);
count_ATGS_low_PA = zeros(1,atg_counter_len);

for i = 1:size(T,1)
        
    orf = char(T{i,"ORF_1"});
    utr5 = char(T{i,"UTR_5"});
    utr5_len_orig = T{i,"UTR5_LEN_ORIG"};
    
    if T{i,"PA1"} >= P_high
        count_ATGS_high_PA = countATGInFrames(utr5, orf, utr5_len_orig, utr_orf_len, count_ATGS_high_PA);
    elseif T{i,"PA1"} <= P_low
        count_ATGS_low_PA = countATGInFrames(utr5, orf, utr5_len_orig, utr_orf_len, count_ATGS_low_PA);
    end
end

figure;
A = [count_ATGS_high_PA(1:3:end); count_ATGS_low_PA(1:3:end)]';

bar(log2(A), 'BarWidth', 1, 'EdgeColor', 'none');
h = get(gca, "Children");
set(h(1), 'FaceColor','g');
set(h(2), 'FaceColor', 'r');

figure;
B = [count_ATGS_high_PA(2:3:end); count_ATGS_low_PA(2:3:end)]';
bar(B, 'BarWidth', 1, 'EdgeColor', 'none');
h = get(gca, "Children");
set(h(1), 'FaceColor','g');
set(h(2), 'FaceColor', 'r');

figure;
C = [count_ATGS_high_PA(3:3:end); count_ATGS_low_PA(3:3:end)]';
bar(C, 'BarWidth', 1, 'EdgeColor', 'none');
h = get(gca, "Children");
set(h(1), 'FaceColor','g');
set(h(2), 'FaceColor', 'r');
