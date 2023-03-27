T = readtable("yeast_parameters_table_with_diff_5utr.xls");


[aa_to_codon, codon_to_aa] = get_aa_to_codon_map();

sample_size = 20;
atg_counter_len = 90;
utr_orf_len = 45;

%use same loop to:
count_ATGS = zeros(1,atg_counter_len);

for i = 1:size(T,1)
    orf = char(T{i,"ORF_1"});
    utr5 = char(T{i,"UTR_5"});
    utr5_len_orig = T{i,"UTR5_LEN_ORIG"};
    count_ATGS = countATGInFrames(utr5, orf, utr5_len_orig, 45, count_ATGS);
end

variablesTypes = repmat({'double'}, 1, atg_counter_len);
randUTRandORFTable = table('Size',[sample_size atg_counter_len], 'VariableTypes', variablesTypes);

count_rands_high_ATG = zeros(1,atg_counter_len);
count_rands_low_ATG = zeros(1,atg_counter_len);

for rand_count = 1:20
    count_ATGS_rand = zeros(1,atg_counter_len);
    for i = 1:size(T,1)
        
        orf = char(T{i,"ORF_1"});
        orf_rand_sample = getRandORF(orf, aa_to_codon);

        utr5 = char(T{i,"UTR_5"});
        
        if length(utr5) > utr_orf_len
            utr5 = utr5(length(utr5)-utr_orf_len+1:end);
        end
        utr5_rand_sample = getRandUTR(utr5);

        utr5_len_orig = T{i,"UTR5_LEN_ORIG"};
        
        count_ATGS_rand = countATGInFrames(utr5_rand_sample, orf_rand_sample,  utr5_len_orig, 45, count_ATGS_rand);
        T.RAND_ORF_D{i} = orf_rand_sample;
    
    end
    count_rands_low_ATG = count_rands_low_ATG + double(count_ATGS_rand <= count_ATGS)
    count_rands_high_ATG = count_rands_high_ATG + double(count_ATGS_rand >= count_ATGS)

    randUTRandORFTable(rand_count, :) = array2table(count_ATGS_rand);
end

count_rands_high_ATG = count_rands_high_ATG ./ 20.0;
count_rands_low_ATG = count_rands_low_ATG ./ 20.0;

position_pMeans = zeros(1,atg_counter_len);
position_pValues = zeros(1,atg_counter_len);

%bigger_than = sum(randUTRandORFTable >= count_ATGS)
%position_pValues(i) = bigger_than / atg_counter_len
for i = 1:atg_counter_len
% 
     bigger_than = sum(randUTRandORFTable{:,i} >= count_ATGS(i));
     smaller_than = sum(randUTRandORFTable{:,i} <= count_ATGS(i));
     position_pValues(i) = min(bigger_than / atg_counter_len, smaller_than / atg_counter_len) ;
%     
     position_pMeans(i) = mean(randUTRandORFTable{:,i});
%     position_std = std(randUTRandORFTable{:,i})
%     position_z = (count_ATGS(i) - position_mean) / position_std
%     position_pValue = 2*normcdf(-abs(position_z),0,1)
% 
%     position_pMeans(i) = position_mean;
%     position_pValues(i) = position_pValue;

end

x = (1:1:32);
y = count_rands_low_ATG < 0.05;
idx = ~isnan(y);
y = y(idx);
y = double(y);
y = 10 * y;

y_green = count_rands_high_ATG < 0.05;
idx = ~isnan(y);
y = y(idx);
y_green = double(y_green);
y_green = 105 * y_green;
idx = ~isnan(y_green);
y_green = y_green(idx);

figure;
A = [count_ATGS(1:3:end); position_pMeans(1:3:end)]';
barColorMap = [0 0 0; 1 1 0];  
%barColorMap(2,:) = [.9 .9 .14];  % Yellow Color for segment 2.
colormap(barColorMap);

bar(log2(A), 'BarWidth', 0.8);

hold on;

y_skip = y(1:3:end);

plot(x(y_skip~=0),y_skip(y_skip~=0),'o','MarkerFaceColor','red','MarkerSize',5);
%plot(x,y_green(1:3:end),'o','MarkerFaceColor','green','MarkerSize',5);

figure;
B = [count_ATGS(2:3:end); position_pMeans(2:3:end)]';
bar(B, 'BarWidth', 0.8);
hold on;
y2 = 10.5 * y(2:3:end);
plot(x(y2~=0),y2(y2~=0),'o','MarkerFaceColor','red','MarkerSize',5);
y2 = y_green(2:3:end);
plot(x(y2~=0),y2(y2~=0),'o','MarkerFaceColor','green','MarkerSize',5);

figure;
C = [count_ATGS(3:3:end); position_pMeans(3:3:end)]';
bar(C, 'BarWidth', 0.8);
hold on;
y3 = 5.5 * y(3:3:end);
plot(x(y3~=0),y3(y3~=0),'o','MarkerFaceColor','red','MarkerSize',5);
y3 = y_green(3:3:end) ./ 2;
plot(x(y3~=0),y3(y3~=0),'o','MarkerFaceColor','green','MarkerSize',5);

%writetable(T, 'yeast_parameters_table_with_diff_5utr.xls')

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



