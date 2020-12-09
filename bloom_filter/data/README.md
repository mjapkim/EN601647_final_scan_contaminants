# Data used for the analysis of bloom filter index for contamination detection
## Mycoplasma 

### Bloom filter index created using this reference genome 
Mycoplasma pneumoniae M129 (mycoplasmas) complete genome:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000027345.1/

### Short read data used for contamination testing
Mycoplasma penumoniae reference strains short read sequencing:
https://www.ncbi.nlm.nih.gov/sra?term=SRP081446

## Synthetic Data
#### Disclaimer: This approach had very low accuracy results for contamination detection with this small genome.

### Bloom filter index creation
Half of the results were created using the contamination reference genome to make the bloom filter index and the other half was made using the target reference genome.

#### Contamination reference genome
contamination_mat_genome.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta

#### Target reference genome
target_mat_genome.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta

### Contamination detection tested on short simulated reads 
simulation.gs100000.cov30.het0.5.err1.contamination4.readsize150.fasta

