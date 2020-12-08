from bloomfilter import BloomFilter
from random import shuffle
import time
import sys


def build_bf(n,p,ref_fasta):
    bloomf = BloomFilter(n,p)
    print("Size of bit array:{}".format(bloomf.size))
    print("False positive Probability:{}".format(bloomf.fp_prob))
    print("Number of hash functions:{}".format(bloomf.hash_count))
    
    mycoplasma_fasta = open(ref_fasta,'r')
    N_count = 0
    read_count = 0
    while True:
        name = mycoplasma_fasta.readline() # read id
        if len(name) == 0:
            break # end of file
        read = mycoplasma_fasta.readline().strip()
        if 'N' not in read:
            bloomf.add(read)
            read_count += 1
        else:
            N_count += 1
    print('N_count = %s'%N_count)
    print('read_count = %s'%read_count)
    mycoplasma_fasta.close()
    return bloomf

def read_seq(sample_name):
    test_fasta = open(sample_name,'r')
    test_words = []
    test_names = []
    N_count = 0
    # read in the test file
    while True:
        name = test_fasta.readline().strip() # name
        if len(name) == 0:
            break # end of file
        read = test_fasta.readline().strip()
        if 'N' not in read:
            test_words.append(read)
            test_names.append(name)
        else:
            N_count += 1
    test_fasta.close()
    print('N_count = %s'%N_count)
    print('test read count = %s'%len(test_words))
    return [test_words, test_names]
    
def detect_contamination_output(file_corr,file_contam,test_words,test_names,bloomf):
    c = open(file_contam,"w")
    m = open(file_corr,"w")
    count = 0
    contam_count = 0
    for i,word in enumerate(test_words):
        if bloomf.check(word):
            count += 1
            m.write(test_names[i])
            m.write('\n')
        else:
            contam_count += 1
            c.write(test_names[i])
            c.write('\n')
    print('Count of exact match words = %s'%count)
    print('Count of contamintant words = %s'%contam_count)
    c.close()
    m.close()

def detect_contamination(test_words,test_names,bloomf):
    count = 0
    contam_count = 0
    for i,word in enumerate(test_words):
        if bloomf.check(word):
            count += 1
        else:
            contam_count += 1
    print('Count of exact match words = %s'%count)
    print('Count of contamintant words = %s'%contam_count)

if __name__=="__main__":
    n = 9080911 #no of items to add (+ 1?) - N?
    fpp = [0.01, 0.05, 0.1, 0.15] #false positive probability
    read_length = 101 # should match both the reference and test sequence
    ref_genome = 'ref_genomes/Mycoplasma_ref.fasta'
    # read in sample to test
    sample_name = '/home/en44704/Desktop/sra_data.fasta'
    outfile_path = '/home/en44704/Desktop/EN601647_final_scan_contaminants/bloom_filter/mycoplas_outputs/'
   
    fpp = [0.15]
    for p in fpp:
        # build bf
        start_build = time.perf_counter()
        bloomf = build_bf(n, p, ref_genome)
        end_build = time.perf_counter()
        print(end_build-start_build)
        # mycoplasma bloom filter index in bloomf

        [test_words, test_names] = read_seq(sample_name)
        #for output files
        file_corr = outfile_path + "bloom_correct_reads" + str(p) + ".txt"
        file_contam = outfile_path + "bloom_contam_reads" + str(p) + ".txt"
        
        # here do the bloomfilter check of membership
        #detect_contamination_output(file_corr,file_contam,test_words,test_names,bloomf)
        
        # without output files for timing
        detect_contamination(test_words,test_names,bloomf)
        
        detected = time.perf_counter()
        print(detected - end_build)
        print("00000000000000000000000000000000000")

