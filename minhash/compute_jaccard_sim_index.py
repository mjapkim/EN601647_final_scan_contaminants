import pickle
import time

len_seed = 10
input1 = open('pkl/simulated_150mers_hash.pkl', 'rb')
genome = pickle.load(input1)

input2 = open('pkl/simulated_short_seq_hash.pkl', 'rb')
short_seq = pickle.load(input2)

# remove reads with same 10 x hash values
# only perform on the 150-mers of genome
# if using more hash functions, should be okay to clean 20,000 short reads
print('orig_g: ', len(genome))
d2 = {tuple(v): k for k, v in genome.items()}
genome = {v: list(k) for k, v in d2.items()}
print('clean_g: ', len(genome))

t0 = time.time()
match_res = {}
# computer Jaccard similarity index for all pairs from set A and set B
# if similarity index = 1, move on to next read comparison
for name_g, hashes_g in genome.items():
    for name_s, hashes_s in short_seq.items():
        tot = len([value for value in hashes_g if value in set(hashes_s)])
        if tot > 0:
            jaccard = float(tot)/len_seed
            if name_s not in match_res:
                match_res[name_s] = jaccard
                if jaccard == 1:
                    # print("match")
                    break
            else:
                if match_res[name_s] < jaccard:
                    match_res[name_s] = jaccard
t1 = time.time()

with open('pkl/simulated_minhash_classify_contam.pkl', 'wb') as output:
    pickle.dump(match_res, output, pickle.HIGHEST_PROTOCOL)

total = t1 - t0
print("classify time: ", total)