from CG_FMIndex import *
import pickle
import time

reads = open('simulation.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta', 'r')

with open('pkl/simulated_input_clean.pkl', 'wb') as output:
    t0 = time.time()
    short_seq_pattern = {}
    N_count = 0
    while True:
        id = reads.readline().strip() # name
        id = id[1:]
        if len(id) == 0:
            break  # end of file
        seq = reads.readline().strip()
        if 'N' not in seq:
            short_seq_pattern[id] = seq
        else:
            N_count += 1
    t1 = time.time()
    pickle.dump(short_seq_pattern, output, pickle.HIGHEST_PROTOCOL)


# print(len(short_seq_pattern))
# print(N_count)

# p = "CGACGAAATTAATACCATCAGGGTATTAAGATGCTACC"
short_seq_matches = {}
with open('pkl/simulated_fm_index.pkl', 'rb') as input:
    t2 = time.time()   
    fm = pickle.load(input)
    for name, seq in short_seq_pattern.items():
        # print("id: ", name)
        matches = sorted(fm.occurrences(seq))
        #print(matches)
        if len(matches) > 0 :
            short_seq_matches[name] = matches
    t3 = time.time()
    
with open('pkl/simulated_pattern_match.pkl', 'wb') as output:
    pickle.dump(short_seq_matches, output, pickle.HIGHEST_PROTOCOL)
del short_seq_matches

total = t1-t0
print("clean up total time: ", total)

total1 = t3-t2
print("classify time: ", total1)
