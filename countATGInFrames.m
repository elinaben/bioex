
function count_ATGS = countATGInFrames(utr5, orf, utr5_len_orig, max_length, count_ATGS)

    count_ATGs_UTR5 = countATGsUTR5(utr5, utr5_len_orig, max_length);
    count_ATGs_ORF = countATGInORF(orf, max_length);
    count_ATGS_new = [count_ATGs_UTR5 count_ATGs_ORF];
    
    count_ATGS = count_ATGS + count_ATGS_new;
end

function count_ATGS = countATGsUTR5(utr5, utr5_len_orig, max_length)
        
    count_ATGS = zeros(1, max_length);
    if isnan(utr5_len_orig) 
        return;
    end
    offset_len = 0;
    if length(utr5) < max_length
        offset_len = max_length - length(utr5);
    elseif length(utr5) > max_length
        %concatenate utr5 to take only the last max_length chars
        utr5 = utr5(length(utr5)-max_length+1:end);
    end 
    
    k = strfind(utr5,"ATG");
    for j=1:numel(k)
        atg_loc = k(j) + offset_len;
        count_ATGS(atg_loc) = count_ATGS(atg_loc) + 1;
    end
end


function count_ATGS = countATGInORF(orf, max_length)
    count_ATGS = zeros(1, max_length);    
    orf = orf(1:min(end,max_length));
    k = strfind(orf,"ATG");
    count_ATGS(k) = count_ATGS(k) + 1;    
end

