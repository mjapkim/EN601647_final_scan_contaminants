from CG_FMIndex import *
import pickle

mycoplasma_fasta = open('sra_data.fasta', 'r')

short_seq_pattern = {}
N_count = 0
while True:
    id = mycoplasma_fasta.readline().strip() # name
    id = id[1:-11]
    if len(id) == 0:
        break  # end of file
    seq = mycoplasma_fasta.readline().strip()
    if 'N' not in seq:
        short_seq_pattern[id] = seq
    else:
        N_count += 1
print(len(short_seq_pattern))
print(N_count)

# p = "CGACGAAATTAATACCATCAGGGTATTAAGATGCTACC"
short_seq_matches = {}
with open('fmIndexed.pkl', 'rb') as input:
    fm = pickle.load(input)
    for name, seq in short_seq_pattern.items():
        # print("id: ", name)
        matches = sorted(fm.occurrences(seq))
        print(matches)
        if len(matches) > 0 :
            short_seq_matches[name] = matches
with open('shortSeq_match.pkl', 'wb') as output:
    pickle.dump(short_seq_matches, output, pickle.HIGHEST_PROTOCOL)
