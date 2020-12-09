import pickle

short_seq = open('/Users/aprilkim/Downloads/simulation.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta', 'r')

def getKmers(seq, size):
    return [kmers["Base"+str(x)].append(seq[x:x+size].upper())for x in range(len(seq) - size + 1)]

def make_kmer_table(seqs, k):
    """ Given read dictionary and integer k, return a dictionary that
    maps each k-mer to the set of names of reads containing the k-mer. """
    table = {}
    # for each 101-mer in genome
    for name, seq in seqs.items():
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in table:
                table[kmer] = set() 
            table[kmer].add(name)
    for key in list(table):
        val = table[key]
        if len(val) < 2:
            del table[key]
    return table

sequence = {}
while True:
    first_line = short_seq.readline() # name
    if len(first_line) == 0:
        break  # end of file
    seq = short_seq.readline().rstrip()
    sequence[first_line] = seq

with open('simulated_pkl/simulated_short_seq_dict.pkl', 'wb') as output:
    pickle.dump(sequence, output, pickle.HIGHEST_PROTOCOL)
del sequence