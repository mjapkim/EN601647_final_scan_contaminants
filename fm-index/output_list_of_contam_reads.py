import pickle

input1 = open('pkl/simulated_pattern_match.pkl', 'rb')
match_tbl = pickle.load(input1)

for name, seq in match_tbl.items():
    print(name)