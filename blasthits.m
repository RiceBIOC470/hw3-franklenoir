function result = blasthits(ascennum,N)
    genbankinfo = getgenbank(ascennum);
    res = {};
    seq = genbankinfo.Sequence;
    [requestID, requestTime] = blastncbi(seq,'blastn','Database','refseq_rna');
    blast_data = getblast(requestID,'WaitTime',requestTime);
    temp = blast_data.Hits(1:N);
    for i = 1:N 
        tempnam = temp(i).Name;
        ref = strsplit(tempnam ,'|');
        res(1,i) = ref(4);
    end
    result = res;
end
