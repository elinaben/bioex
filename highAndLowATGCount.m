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

    if length(utr5) > utr_orf_len
        utr5 = utr5(length(utr5)-utr_orf_len+1:end);
    end
    if T{i,"PA1"} >= P_high
        count_ATGS_high_PA = countATGInString(utr5, orf, count_ATGS_high_PA);
    elseif T{i,"PA1"} <= P_low
        count_ATGS_low_PA = countATGInString(utr5, orf, count_ATGS_low_PA);
    end
end

figure;
A = [count_ATGS_high_PA(1:3:end); count_ATGS_low_PA(1:3:end)]';
barColorMap = [0 0 0; 1 1 0];  
%barColorMap(2,:) = [.9 .9 .14];  % Yellow Color for segment 2.
colormap(barColorMap);

bar(log2(A), 'BarWidth', 0.8);

figure;
B = [count_ATGS_high_PA(2:3:end); count_ATGS_low_PA(2:3:end)]';
bar(B, 'BarWidth', 0.8);

figure;
C = [count_ATGS_high_PA(3:3:end); count_ATGS_low_PA(3:3:end)]';
bar(C, 'BarWidth', 0.8);

function count_ATGS = countATGInString(utr5, orf, count_ATGS)

    utr_orf_len = 90;
    offset_len = 0;
    if length(utr5) < utr_orf_len
        offset_len = utr_orf_len - length(utr5);
    end

    utr5_and_orf = strcat(utr5, orf(1:min(end,utr_orf_len)));
    k = strfind(utr5_and_orf,"ATG");
    for j=1:numel(k)
        atg_loc = k(j) + offset_len;
        if atg_loc > 180
            disp("Loc 180" + utr5);
            disp(orf);
            disp(atg_loc);
            disp(k(j));
            continue;
        end
        count_ATGS(atg_loc) = count_ATGS(atg_loc) + 1;
    end
end