# EN601647_final_scan_contaminants

<b>Summary Goal:</b> Identify Mycoplasma genome contamination in the NCBI short read dataset for S.Cerevisiae through approximate matching. Compare three different indexing methods to index Mycoplasma genome as contaminants, evaluating them in terms of speed and accuracy of contamination location identification, and report performance.

<b>What is the method you want to develop?</b>
High-throughput sequencing technologies have been used to create genomes of different species in a cost-effective and rapid manner. Sequences obtained, however, may be contaminated with DNA from sources other than the sample and different methods to identify and remove contamination from genomic datasets have been developed. Current techniques use an offline approach, where the given sample is compared to target datasets to identify matches.

<b>Clarification of ‘kmer profiling’:</b>
In our project, we would like to approach this problem by evaluating the k-mer profile for the contaminant reference and compare it to different samples. By k-mer profile, we mean comparing the coverage and frequency of k-mers that appear in the reference to the sample. We will use the expectation that a contaminated sample would have slightly less defined heterozygotic peaks. K-mers that do not fit in the reference profile or within the sample error peaks could indicate contamination. We can identify k-mers that appear in our contaminant reference as an indication for contamination.

We will compare this k-mer profiling approach against three other indexing techniques: Bloom filters, FM index, and minhash.
