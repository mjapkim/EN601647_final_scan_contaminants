from bloomfilter import BloomFilter
from random import shuffle

n = 9000000 #no of items to add (+ 1?) - N?
p = 0.05 #false positive probability

bloomf = BloomFilter(n,p)
print("Size of bit array:{}".format(bloomf.size))
print("False positive Probability:{}".format(bloomf.fp_prob))
print("Number of hash functions:{}".format(bloomf.hash_count))

#mycoplasma_fasta = open('portion.fasta','r')
mycoplasma_fasta = open('/home/en44704/Desktop/sra_data.fasta','r')

reads_present = []
N_count = 0
while True:
    name = mycoplasma_fasta.readline() # name
    if len(name) == 0:
        break # end of file
    read = mycoplasma_fasta.readline().strip()
    if 'N' not in read:
        bloomf.add(read)
    else:
        N_count += 1

print('N_count = %s'%N_count)
# mycoplasma bloom filter index
