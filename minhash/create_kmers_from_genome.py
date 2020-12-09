import pickle

fasta = open('/Users/aprilkim/Downloads/target_mat_genome.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta', 'r')

def getKmers(seq, size):
    return [[seq[x:x+size].upper(), "base"+str(x)] for x in range(len(seq) - size + 1)]

genome = []
next(fasta) # header starting with <
while True:
    seq = fasta.readline().strip()
    if len(seq) == 0:
        break  # end of file
    genome.append(seq)
genome = ''.join([x for x in genome])

kmers = getKmers(genome, 150)
kmers_dict = {}
for line in kmers:
    kmers_dict[line[1]] = line[0]

with open('simulated_150mers_dict.pkl', 'wb') as output:
    pickle.dump(kmers_dict, output, pickle.HIGHEST_PROTOCOL)
del kmers_dict