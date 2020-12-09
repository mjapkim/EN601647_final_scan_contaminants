from CG_FMIndex import *
import pickle
import time

fasta = open("target_mat_genome.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta", "r")

genome = []
next(fasta) # header starting with <
while True:
    seq = fasta.readline().strip()
    if len(seq) == 0:
        break  # end of file
    genome.append(seq)
genome = ''.join([x for x in genome])

with open('pkl/simulated_fm_index.pkl', 'wb') as output:
    t0 = time.time()
    fm = FmIndex(genome)
    t1 = time.time()
    total = t1 - t0
    pickle.dump(fm, output, pickle.HIGHEST_PROTOCOL)

del fm

print(total)