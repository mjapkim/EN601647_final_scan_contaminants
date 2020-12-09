import pickle
import time

fasta = open('target_mat_genome.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta', 'r')
short_seq = open('simulation.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta', 'r')

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

t0 = time.time()
kmers = getKmers(genome, 150)
kmers_dict = {}
for line in kmers:
    kmers_dict[line[1]] = line[0]
t1 = time.time()

with open('pkl/simulated_150mers_dict.pkl', 'wb') as output:
    pickle.dump(kmers_dict, output, pickle.HIGHEST_PROTOCOL)
del kmers_dict

total = t1 - t0
print("time taken to create 150-mers: ", total)

t2 = time.time()
sequence = {}
while True:
    first_line = short_seq.readline() # name
    if len(first_line) == 0:
        break  # end of file
    seq = short_seq.readline().rstrip()
    sequence[first_line] = seq
t3 = time.time()

with open('pkl/simulated_short_read_dict.pkl', 'wb') as output:
    pickle.dump(sequence, output, pickle.HIGHEST_PROTOCOL)
del sequence

total1 = t3 - t2
print("time taken to clean short reads: ", total1)