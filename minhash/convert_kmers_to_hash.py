import pickle
import mmh3
import time

def string_to_min_hash(Astr, k, seed=0):
    return min([mmh3.hash(Astr[i:i+k], seed) for i in range(len(Astr)-k+1)])

# def jaccard_min_kmer_hash(Astr, Bstr, k):
#     return 1 if string_to_min_hash(Astr, k) == string_to_min_hash(Bstr, k) else 0

input0 = open('pkl/simulated_150mers_dict.pkl', 'rb')
sequence = pickle.load(input0)

input1 = open('pkl/simulated_short_read_dict.pkl', 'rb')
sequence1 = pickle.load(input1)

t0 = time.time()
hash_dict_genome = {}
for name, seq in sequence.items():
    str_hash = [string_to_min_hash(seq, 60, k) for k in range(10)]
    hash_dict_genome[name] = str_hash
t1 = time.time()

t2 = time.time()
hash_dict_short = {}
for name, sequence in sequence1.items():
    str_hash = [string_to_min_hash(seq, 60, k) for k in range(10)]
    hash_dict_short[name] = str_hash
t3 = time.time()

with open('pkl/simulated_150mers_hash.pkl', 'wb') as output:
    pickle.dump(hash_dict_genome, output, pickle.HIGHEST_PROTOCOL)

with open('pkl/simulated_short_seq_hash.pkl', 'wb') as output:
    pickle.dump(hash_dict_short, output, pickle.HIGHEST_PROTOCOL)

total = t1 - t0
print("time taken convert 150-mers to hash: ", total)
total1 = t3 - t2
print("time taken to convert reads to hash: ", total)