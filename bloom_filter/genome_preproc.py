import time
import sys

def remove_newline(s):
    char = '\n'
    count_nl = s.count(char)
    s = list(s)
    while count_nl:
        s.remove(char)
        count_nl -= 1
    return ''.join(s)

def circular_read(a, k, ind):
    i = ind
    n = len(a)
    read = ""
    while i < k + ind: # when we have gone + k chars forward
        read = read + a[(i % n)]
        i += 1
    return read

if __name__=="__main__":
    read_length = 150

    tic = time.perf_counter()
    
    fasta = open('/home/en44704/Desktop/mycoplasma_pneumoniae_complete_genome.fasta','r')
    #fasta = open('/home/en44704/Desktop/contamination_mat_genome.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta','r')
    #fasta = open('/home/en44704/Desktop/target_mat_genome.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta','r')
    
    text_file = open("Mycoplasma_ref.fasta","w")
    #text_file = open("synthetic_contam_ref.fasta","w")
    #text_file = open("synthetic_target_ref.fasta","w")
    read_count = 0
    N_count = 0
    reads = ""
    header = fasta.readline() # header
    fake_reads = fasta.read()
    # take out all the \n 
    reads = remove_newline(fake_reads)
    rem_new = time.perf_counter()
    print(rem_new-tic)
    for i,char in enumerate(reads): 
        read = circular_read(reads,read_length,i)
        name = str(i) + 'read'
        if 'N' not in read:
            text_file.write(name)
            text_file.write('\n')
            text_file.write(read)
            text_file.write('\n')
            read_count +=1
        else:
            N_count += 1
    toc = time.perf_counter()
    text_file.close()
    
    print('read_count = %s'%read_count)
    print('N_count = %s'%N_count)
    print(toc-tic)
    fasta.close()
