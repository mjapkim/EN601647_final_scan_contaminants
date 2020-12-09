import pickle5 as pickle
import time

input1 = open('simulated_pkl/simulated_minhash_classify_contam.pkl', 'rb')
res = pickle.load(input1)

clean = {}
contam = {}
print(len(res))
for key, jaccard_sim in res.items():
    if jaccard_sim != 1:
        contam[key] = jaccard_sim
    else:
        clean[key] = jaccard_sim
print(len(contam))
print(len(clean))
