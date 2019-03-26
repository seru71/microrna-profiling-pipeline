
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq

#
# extract mirna
#

def extract_mirna_sequences(input_file_handle, output, prefix='hsa-', to_dna=True):
    for r in SeqIO.parse(input_file_handle, "fasta"):
        if not prefix or (prefix and r.id.startswith(prefix)):
            if to_dna:
                r.seq = r.seq.back_transcribe()
            SeqIO.write(r, output, "fasta")   
                #output.write(">" + r.description + " cDNA\n")
                #output.write(str(r.seq.back_transcribe()) + '\n')
            #else:
                #output.write(">" + r.description)
                #output.write(str(r.seq) + '\n')


#
# collapsing
#

def collapse_identicals(input_file_handle, output, collapse_rc_identicals=False):
    """ join records with identical sequences, optionally considering reverse complements """
    mirnas = {}
    seqid=""
    for r in SeqIO.parse(input_file_handle, "fasta"):
        if str(r.seq) in mirnas: 
            mirnas[str(r.seq)].append(r)
        elif collapse_rc_identicals and str(r.seq.reverse_complement()) in mirnas:
            r_rc = r.reverse_complement()
            r_rc.id = r.id #+"_RC"
            r_rc.description = r.description + "reverse complement"
            mirnas[str(r_rc.seq)].append(r_rc)
        else: mirnas[str(r.seq)] = [r]

    for rs in mirnas.values():
       if len(rs) > 1:
           ids = [r.id for r in rs]
           mimats = [r.description.split(' ')[1] for r in rs]
           rs[0].id = "_".join(ids)
           rs[0].description = "_".join(mimats)
           
       SeqIO.write(rs[0], output, "fasta")


#
# padding
#

def pad_sequence(s, left_n, right_n, char):
    return s.rjust(len(s)+left_n, char).\
             ljust(len(s)+left_n+right_n, char)

def pad_fasta_records(input_file_handle, output, left_padding_len, right_padding_len, nucleotide='N'):
    """ Pads fasta records """
    for r in SeqIO.parse(input_file_handle, "fasta"):
        r.seq = Seq(pad_sequence(str(r.seq), left_padding_len, right_padding_len, nucleotide))
        SeqIO.write(r, output, "fasta")
    
    
#
# fasta2fastq
#

def mirbase_fasta_to_fastq(fasta_handle, output):
    for r in SeqIO.parse(fasta_handle, "fasta"):
        r.letter_annotations["phred_quality"] = [40] * len(r)
        SeqIO.write(r, output, "fastq")


#
# collapse multimapers
#

def getSAMField(name, split_line):
    FIELDS = {'query':0, 'flag':1, 'ref':2, 'mq':4}
    return split_line[FIELDS[name]]

def getQuery(split_line):
    return getSAMField('query', split_line)

def getFlag(split_line):
    return int(getSAMField('flag', split_line))
    
def getRef(split_line):
    return getSAMField('ref', split_line)
    
def getMQ(split_line):
    return int(getSAMField('mq', split_line))
    
def hasExtraTags(split_line):
    return len(split_line)>11
    
def getTagValue(tag, split_line):
    return [e for e in split_line[11:] if e.startswith(tag)][0][len(tag+":"):]
    
def isHeader(line):
    return line[0] == "@"

def numMismatches(split_line):
    return int(getTagValue('NM:i',split_line))

def getXAs(split_line):
    xa = getTagValue("XA:Z", split_line)
    xa_elems   = [e.split(',')[0] for e in xa.split(';') if e != ""]
    xa_revcomp = [int(e.split(',')[1])<0 for e in xa.split(';') if e != ""]
    # store orientations of the sequeneces in a dict
    return {s:o for s,o in zip(xa_elems, xa_revcomp)}

def get_multimapping_clusters(input_file):
    """ check alignment source and parse accordingly """
    for l in input_file:
        ls = l.strip().split('\t')
        if not isHeader(l):
            raise Exception("No recognized tool found in the SAM header. Supported alignments are bowtie2 and BWA aln")
        if ls[0]=="@PG" and ls[1]=="ID:bwa":
            return get_multimapping_clusters_bwa(input_file)
        elif ls[0]=="@PG" and ls[1]=="ID:bowtie2":
            return get_multimapping_clusters_bowtie(input_file)

def get_multimapping_clusters_bowtie(input_file):
    """ only one-direction mapping is supported. for bidriectional mapping of reads, use BWA """
    direction=False
    clusters = []
    id_cluster_map = {}

    for l in input_file:
        
        ls = l.strip().split('\t')
        if isHeader(l) or getMQ(ls)>1 or numMismatches(ls)>0: # will all multimappers get MQ=1?
            continue
        
        query = getQuery(ls)
        ref = getRef(ls)
                
        already_assigned = query in id_cluster_map, ref in id_cluster_map
        
        if all(already_assigned):
            if id_cluster_map[query] == id_cluster_map[ref]: 
                continue
            
            # both query and ref are members of distinct existing clusters - merge these clusters
            clusters[id_cluster_map[ref][0]] = clusters[id_cluster_map[ref][0]].\
                                                union(clusters[id_cluster_map[query][0]])
            clusters[id_cluster_map[query][0]] = None
            id_cluster_map[query] = (id_cluster_map[ref][0], direction)
        elif not any(already_assigned):
            # none of the elements of s are already in a cluster, make a new one
            id_cluster_map[query] = (len(clusters), direction)
            id_cluster_map[ref] = (len(clusters), direction)
            clusters.append(set([query, ref]))
            #print 'created new cluster', s
        else:
            # either query or ref are members of an existing cluster, add the other
            if query in id_cluster_map:
                clusters[id_cluster_map[query][0]].add(ref)
                id_cluster_map[ref] = (id_cluster_map[query][0], direction)
            else:
                clusters[id_cluster_map[ref][0]].add(query)
                id_cluster_map[query] = (id_cluster_map[ref][0], direction)
    
    # remove Nones and keep only orientations in the dict
    return [c for c in clusters if c is not None], \
            {e: id_cluster_map[e][1] for e in id_cluster_map}



def get_multimapping_clusters_bwa(input_file):
    clusters = []
    id_cluster_map = {}

    for l in input_file:
        ls = l.strip().split('\t')
        if isHeader(l) or getMQ(ls)>0 or not hasExtraTags(ls): 
            continue
        
        # extract query ID and query orientation
        query = getQuery(ls)
        query_revcomp = bool((getFlag(ls) >> 4) & 1)
        
        # ref id
        ref = getRef(ls)
        
        # extract alternative hits and their orientation to a dict
        s_revcomp = getXAs(ls)
        #s_revcomp[query] = False 
        s_revcomp[ref] = query_revcomp # revcomp reference if query mapped reverse
        
        # make a set of similar seqeuences (to skip duplicated ids)
        s = set([ref] + s_revcomp.keys())
        
        # check which sequences have been already assigned to a cluster
        already_assigned = [e in id_cluster_map for e in s]
        
        if all(already_assigned):   # skip as this is query,ref & ref,query situation
            continue
        elif not any(already_assigned):
            # none of the elements of s are already in a cluster, make a new one
            index_of_the_new_cluster = len(clusters)
            for e in s:
                id_cluster_map[e] = (index_of_the_new_cluster, s_revcomp[e])
            clusters.append(s)
            #print 'created new cluster', s
        else:
            # some of the elements of s are members of an existing cluster - merge s with these clusters
            
            cluster_ids = set([id_cluster_map[e][0] for e in s if e in id_cluster_map])
            
            cluster_id = cluster_ids.pop()                     
            #print 'adding', s_revcomp, '\nTo cluster', [(sid,id_cluster_map[sid][1]) for sid in clusters[cluster_id]], 'because', already_assigned
            
            # get all sequence_ids in s that have already beed assigned this cluster
            overlap_ids = s.intersection(clusters[cluster_id])
            seq_id = next(iter(overlap_ids))
            reverse_orientations = False
            if id_cluster_map[seq_id][1] is s_revcomp[seq_id]:
                # if orientations of the overlapping sequence in the existing cluster
                # and the added set are identical, we don't modify orientation flags, and skip adding
                #print 'orientations match', [(sid,id_cluster_map[sid][1]) for sid in s.intersection(clusters[cluster_id])], 'vs', s_revcomp
                for e in overlap_ids:
                    if id_cluster_map[e][1] is not s_revcomp[e]:
                        raise Exception("Orientations don't match")
                    #id_cluster_map[e] = (cluster_id, s_revcomp[e])
            else:
                # if orientations of the overlapping sequence in the existing cluster
                # and the added set are opposite, we include newly added sequences with reversed orientation flag
                #print 'orientations dont match', [(sid,id_cluster_map[sid][1]) for sid in s.intersection(clusters[cluster_id])], 'vs', s_revcomp
                for e in overlap_ids:
                    if id_cluster_map[e][1] is s_revcomp[e]:
                        raise Exception("Orientations don't match")
                    #id_cluster_map[e] = (cluster_id, not s_revcomp[e])
                reverse_orientations = True
            
            s -= overlap_ids
            
            # other clusters (if there were any)
            for other_cluster_id in cluster_ids:
                
                overlap_ids = s.intersection(clusters[other_cluster_id])
                seq_id = next(iter(overlap_ids))
                
                if (id_cluster_map[seq_id][1] is s_revcomp[seq_id] and not reverse_orientations) or \
                   (id_cluster_map[seq_id][1] is not s_revcomp[seq_id] and reverse_orientations):
                    for k,v in id_cluster_map.items():
                        if v[0] == other_cluster_id:
                            id_cluster_map[k] = (cluster_id, id_cluster_map[k][1])
                else:
                    for k,v in id_cluster_map.items():
                        if v[0] == other_cluster_id:
                            id_cluster_map[k] = (cluster_id, not id_cluster_map[k][1])
                                
                # absorb the other cluster
                clusters[cluster_id] = clusters[cluster_id].union(clusters[other_cluster_id])
                clusters[other_cluster_id] = None
                
                s -= overlap_ids
                            
            for e in s:
                id_cluster_map[e] = (cluster_id, not s_revcomp[e] if reverse_orientations else s_revcomp[e])
                clusters[cluster_id].add(e)
                
            # merge clusters
            clusters[cluster_id] = clusters[cluster_id].union(s)

    # remove Nones and keep only orientations in the dict
    return [c for c in clusters if c is not None], \
            {e: id_cluster_map[e][1] for e in id_cluster_map}
    


def choose_seed(sequences, start_seed_len=10):
    """ take inner sequence occuring only once in every sequence """
    lengths = [len(s) for s in sequences]
    shortest = sequences[lengths.index(min(lengths))]
    
    for seed_len in range(min(start_seed_len, len(shortest)), len(shortest)+1):
        start = (len(shortest)-seed_len)/2
        seed = shortest[start:(start+seed_len)]
        # check if all seqeuences have exactly one occurence of seed        
        if all([len(re.findall(seed, s)) == 1 for s in sequences]): 
            return seed
    raise Exception('No unique seed could be determined for %s' % sequences)
    
def align_sequences(sequences):
    """ aligns sequences by shifting """
    seed = choose_seed(sequences)
    seed_starts = [seq.find(seed) for seq in sequences]
    innerst = max(seed_starts)
    shifted = [seq.rjust(len(seq) + innerst - seed_start, " ") for seq, seed_start in zip(sequences, seed_starts)]
    max_length = max([len(s) for s in shifted])
    return [s.ljust(max_length, " ") for s in shifted]
    
def get_consensus(sequences):
    """ for aligned sequences get a consensus sequence """
    
    # all seqeunces are equal in length
    assert len(set([len(s) for s in sequences])) == 1
    
    def consensus_char(chars):
        s = set(chars)
        if " " in s: 
            s.remove(" ")
        return s.pop() if len(s)==1 else 'N'
    
    consensus = ""
    for i in range(0, len(sequences[0])):
        nucs = [s[i] for s in sequences]
        consensus += consensus_char(nucs)
    
    return consensus
    
    
def cluster_multimappers(mapping_file_handle, output, reference_fasta_file, **other_args):
    """ join records with overlapping sequence """
    
    def print_fasta(out, seq_id, seq):
        out.write(">" + seq_id + "\n" + seq + "\n")
    
    from pyfaidx import Fasta
    ref = Fasta(reference_fasta_file)
    remaining_sequences = ref.keys()

    clusters, revcomp_flags = get_multimapping_clusters(mapping_file_handle, **other_args) 
    for cluster in clusters:
        seqs = []
        for seq_id in cluster:
            s = ref[seq_id][:]
            if revcomp_flags[seq_id]:
                #s.name = s.name+" (revcomped)"
                s.seq = s.reverse.complement.seq
            seqs.append(s)
            remaining_sequences.remove(seq_id)
        
        new_id = "_".join(cluster)+"_consensus"
        cluster_consensus = get_consensus(align_sequences([s.seq for s in seqs]))
        print_fasta(output, new_id, cluster_consensus)

    for seq_id in remaining_sequences:
        print_fasta(output, seq_id, ref[seq_id][:].seq)

    
def print_clustered(mapping_file_handle, output, reference_fasta_file, **other_args):
    """ reference fasta must be indexed """
    
    clusters, revcomp_flags = get_multimapping_clusters(mapping_file_handle, **other_args)
    
    from pyfaidx import Fasta
    ref = Fasta(reference_fasta_file)

    for cluster in clusters:
        seqs = []
        for seq_id in cluster:
            s = ref[seq_id][:]
            if revcomp_flags[seq_id]:
                s.name = s.name+" (revcomped)"
                s.seq = s.reverse.complement.seq
            seqs.append(s)
        
        aligned_seqs = align_sequences([s.seq for s in seqs])
        for i,s in enumerate(aligned_seqs):
            output.write(s + " " + seqs[i].name + "\n")
        output.write("-----------\n")
        output.write(get_consensus(aligned_seqs) + "\n")
        output.write("+++++++++++\n\n")



#
# ------------------------------------------------------------------#
#


def perform_task(task, input_file_handle, args):
    if task == 'extract-hsa':
        extract_mirna_sequences(input_file_handle, sys.stdout)
    elif task == 'collapse-identical':
        collapse_identicals(input_file_handle, sys.stdout)
    elif task == 'pad':
        pad_fasta_records(input_file_handle, sys.stdout, \
                          args.left_padding, args.right_padding, \
                          args.padding_char)
    elif task == 'fasta2fastq':
        mirbase_fasta_to_fastq(input_file_handle, sys.stdout)
    elif task == 'collapse-multimapers':
        #print_clustered(input_file_handle, sys.stdout, str(args.reference_fasta))
        cluster_multimappers(input_file_handle, sys.stdout, str(args.reference_fasta))

"""
Workflow:

wget mature.fa.gz

zcat mature.fa.gz \
 | python mirbase_prep.py -t extract-hsa \
 | python mirbase_prep.py -t collapse-identical > mature.hsa.collapsed.fa \
&& samtools faidx mature.hsa.collapsed.fa \
&& python mirbase_prep.py -t fasta2fastq -i mature.hsa.collapsed.fa > mature.hsa.collapsed.fq \ 
&& bowtie2-build mature.hsa.collapsed.fa mature.hsa.collapsed.fa

bowtie2 --norc -k100 -L11 -x mature.hsa.collapsed.fa mature.hsa.collapsed.fq 2>/dev/null \
 | python mirbase_prep.py -t collapse-multimapers -r mature.hsa.collapsed.fa > mature.hsa.multcollapsed.fa
 
bowtie2-build mature.hsa.multcollapsed.fa mature.hsa.multcollapsed.fa

... and mature.hsa.multcollapsed.fa can be used for mapping.

"""

if __name__ == '__main__':
    
    from optparse import OptionParser
    parser = OptionParser(version="%prog 1.0", usage = "\n\n %prog [extract-hsa|collapse-identical|pad|fasta2fastq|collapse-multimapers] input_file")
    parser.add_option("--task",            "-t", type='string', action="store", dest="task",            help="Task to perform")
    parser.add_option("--input_file",      "-i", type='string', action="store", dest="input_file",      help="Input file for the task [stdin]")
    parser.add_option("--reference_fasta", "-r", type='string', action="store", dest="reference_fasta", help="Reference file for the task. Required for tasks: cluster")
    parser.add_option("--left_padding",    "--pl",      type='int',    action="store", dest="left_padding",    default=3,   help="Padding length on the left. Required for tasks: pad")
    parser.add_option("--right_padding",   "--pr",      type='int',    action="store", dest="right_padding",   default=3,   help="Padding length on the right. Required for tasks: pad")
    parser.add_option("--padding_nucleotide", "--pn",   type='string', action="store", dest="padding_char",    default='N', help="Padding nucleotide character [N]. Required for tasks: pad")
    (args, _) = parser.parse_args()
    
    implemented_tasks = ['extract-hsa','collapse-identical', 'pad', 'fasta2fastq', 'collapse-multimapers']
    if args.task not in implemented_tasks:
        parser.print_help()
        sys.exit(1)

    stdin_aliases = ["-", "stdin", "/dev/stdin"]
    if args.input_file in stdin_aliases or args.input_file is None:
        perform_task(args.task, sys.stdin, args)
    else:
        with open(args.input_file) as input_file: 
            perform_task(args.task, input_file, args)

    

