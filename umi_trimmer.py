"""
UMI trimming module.
Pawel Sztromwasser, 04.2018
"""


"""
function by lh3 and brenp, downloaded from https://github.com/lh3/readfq/blob/master/readfq.py
"""
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


import tre

def find_adapter(read, adapter, error):
    fz = tre.Fuzzyness(maxerr=error)
    pt = tre.compile(adapter)
    matches = pt.search(read, fz)
    if matches is None:
        return None
    else:  # return first match
        return(matches.groups()[0])
    

def openfile(filename, mode):
    import gzip
    return gzip.open(filename, mode) if filename[-3:]=='.gz' else open(filename, mode)


def trim_umi(input_fname, output_fname, preumi_adapter, 
		adapter_mismatch=2, min_umi_len=12, max_umi_len=12, min_read_len=15):
    
    cnt_ok, cnt_short, cnt_no_umi, cnt_short_umi = 0, 0, 0 ,0  
    with openfile(input_fname, 'rU') as infq, openfile(output_fname, 'w') as outfq:
        
        for name, seq, qual in readfq(infq):
            
            match = find_adapter(seq, preumi_adapter, adapter_mismatch)
        
            if match==None:  # drop the read - no adapter = no UMI
                cnt_no_umi+=1
                continue
        
            start, end = match 
            
            if start < min_read_len: # read too short - drop it
                cnt_short+=1
                #print "dropping short read %s nt" % str(start)
                continue
                
            if len(seq) - end < min_umi_len: # UMI shorter than minimal
                cnt_short_umi+=1
               #print "dropping read with UMI len %s" % str(len(seq) - end)
                continue
            
            cnt_ok+=1
            umi_len = min(max_umi_len, len(seq) - end)
            umi = seq[end:(end+umi_len)]
            name += "_UMI:"+umi
            seq = seq[:start]
            qual = qual[:start]
            
            outfq.write('\n'.join(["@"+name,seq,"+",qual]) + '\n')
    
    return {'OK_READS':cnt_ok, 
            'SHORT_READS': cnt_short,
            'NO_UMI': cnt_no_umi,
	    'SHORT_UMI': cnt_short_umi}
    

