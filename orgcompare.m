function [humanhit,organismhit] = orgcompare(ascennum)
    res = blasthits(ascennum,50);
    humnum = 51;
    orgnum = 51;
    humanseq = '';
    orgseq = '';
    originorg = getgenbank(char(res(1)));
    origin = strtrim(originorg.SourceOrganism(1,:));
    
    for i = 2:50 % first hit is original searched term
        if orgnum > 50 || humnum > 50
            temp = getgenbank(char(res(i)));
            if contains(temp.Source,'Homo sapiens','IgnoreCase',true) && i < humnum
                humanseq = char(res(i));
                humnum = i;
            elseif ~contains(temp.Source,origin,'IgnoreCase',true) && i < orgnum
                orgseq = char(res(i));
                orgnum = i;
            end
        else
            break
        end
    end
    if orgnum > 50 
        disp("No other organism in top 50 hits");
    end
    if humnum > 50 
        disp("No other human hit in top 50 hits");
    end
    humanhit = humanseq;
    organismhit = orgseq;
end