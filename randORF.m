T = readtable("yeast_parameters_table_with_diff_5utr.xls");


[aa_to_codon, codon_to_aa] = get_aa_to_codon_map();

sample_size = 20;
atg_counter_len = 96;
utr_orf_len = 45;

%use same loop to:
count_ATGS = zeros(1,atg_counter_len);

for i = 1:size(T,1)
        
    orf = char(T{i,"ORF_1"});
    utr5 = char(T{i,"UTR_5"});
    utr5_len_orig = T{i,"UTR5_LEN_ORIG"};

    if length(utr5) > utr_orf_len
        utr5 = utr5(length(utr5)-utr_orf_len+1:end);
    end

    count_ATGS = countATGInString(utr5, orf, count_ATGS);
end

variablesTypes = repmat({'double'}, 1, atg_counter_len);
randUTRandORFTable = table('Size',[sample_size atg_counter_len], 'VariableTypes', variablesTypes);

count_rands_high_ATG = zeros(1,atg_counter_len);
count_rands_low_ATG = zeros(1,atg_counter_len);

for rand_count = 1:sample_size
    count_ATGS_rand = zeros(1,atg_counter_len);
    for i = 1:size(T,1)
        
        orf = char(T{i,"ORF_1"});
        orf_rand_sample = get_random_orf(orf, aa_to_codon);

        utr5 = char(T{i,"UTR_5"});
        
        if length(utr5) > utr_orf_len
            utr5 = utr5(length(utr5)-utr_orf_len+1:end);
        end
        utr5_rand_sample = get_random_utr(utr5);
        
        count_ATGS_rand = countATGInString(utr5_rand_sample, orf_rand_sample, count_ATGS_rand);
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

function rand_orf = get_random_orf(orf, aa_to_codon)
    
    aa_keys = keys(aa_to_codon);
    rand_orf = orf(1:min(90,length(orf)));
    for j = 1:numel(aa_keys)
        codons = aa_to_codon(char(aa_keys(j)));
        idxes = [];
        for k=1:numel(codons)
            codon_idxes = strfind(orf, codons(k));
            codon_idxes = codon_idxes(mod(codon_idxes,3) == 1);
            idxes = [idxes, codon_idxes];
        end

        %generate a random permutation of the indices
        perm_idxes = randperm(length(idxes));
        idxes_after_perm = idxes(perm_idxes);

        for k=1:numel(idxes_after_perm)
            pre_idx = idxes(k);
            new_idx = idxes_after_perm(k);
            if strcmp(orf(pre_idx:pre_idx+2), orf(new_idx:new_idx+2)) == 0
                if mod(pre_idx-1,3) == 0 % 
                    rand_orf(new_idx:new_idx+2) = orf(pre_idx:pre_idx+2);
                end
            end
        end
    end
end


function ran_utr = get_random_utr(utr)
    n = length(utr);
    %generate a random permutation of the indices
    idx = randperm(n);
    ran_utr = utr(idx);
end

function count_ATGS = countATGInString(utr5, orf, count_ATGS)

    utr_orf_len = 45;
    offset_len = 0;
    if length(utr5) < utr_orf_len
        offset_len = utr_orf_len - length(utr5);
    end

    utr5_and_orf = strcat(utr5, orf(1:min(end,51)));
    k = strfind(utr5_and_orf,"ATG");
    for j=1:numel(k)
        atg_loc = k(j) + offset_len;
        if atg_loc > 95
            disp("Loc 96" + utr5);
            disp(orf);
            disp(atg_loc);
            disp(k(j));
            continue;
        end
        count_ATGS(atg_loc) = count_ATGS(atg_loc) + 1;
    end
end
