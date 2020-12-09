import pickle
import mmh3
import time

def string_to_min_hash(Astr, k, seed=0):
    return min([mmh3.hash(Astr[i:i+k], seed) for i in range(len(Astr)-k+1)])

# def jaccard_min_kmer_hash(Astr, Bstr, k):
#     return 1 if string_to_min_hash(Astr, k) == string_to_min_hash(Bstr, k) else 0

input0 = open('simulated_pkl/simulated_150mers_dict.pkl', 'rb')
sequence = pickle.load(input0)

input1 = open('simulated_pkl/simulated_short_seq_dict.pkl', 'rb')
sequence1 = pickle.load(input1)


hash_dict = {}
for name, seq in sequence.items():
    str_hash = [string_to_min_hash(seq, 60, k) for k in range(10)]
    hash_dict[name] = str_hash

hash_dict1 = {}
for name, sequence in sequence1.items():
    str_hash = [string_to_min_hash(seq, 60, k) for k in range(10)]
    hash_dict1[name] = str_hash

with open('simulated_pkl/simulated_150mers_hash.pkl', 'wb') as output:
    pickle.dump(hash_dict, output, pickle.HIGHEST_PROTOCOL)

with open('simulated_pkl/simulated_short_seq_hash.pkl', 'wb') as output:
    pickle.dump(hash_dict1, output, pickle.HIGHEST_PROTOCOL)