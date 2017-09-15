# Scripts used for analyzing SMRT-cappable-seq data
Developped and maintained by Bo Yan (New England Biolabs). 

All python scripts are based on python 2.7.

See details for each function in README.txt.

- ### pacbio_trim.py
```
Require python regex module
Usage:
python pacbio_trim.py length --input input.fasta --fa
python pacbio_trim.py filer --input --output --filter
python pacbio_trim.py poly --input --output --filter
```
Description:

**length** function analyzes the length of reads in either fastq or fasta file.

**filter** function removes the reads with concatamers and converts all the reads into the same direction as RNA.

**poly** function remove 3'-polyA tail and 5'-polyC adapter sequences in the reads, also extracts UID.

The trimmed reads could be used for mapping using BlazR.


- ### adjust_3end.py
```
Require bedtools
Usage:
python adjust_3end.py --input --output --refer --path
```
Description:

This function adjusts the 3'end for the mapped reads in the bed file (converted from bam file), which removes the mapping due to the addition of polyA tail.

- ### TSS_analysis.py
```
Usage:
python TSS_analysis.py count --input inputbed --output count.output
python TSS_analysis.py cluster --input count.output --output --control --cutoff (default5)
```
Description:

**count** function counts the number of reads starting at the same 5'end and ending at the same 3'end.

**cluster** function cluters the nearby TSSs, whose distance is <=cutoff.

- ### binomialtest.R
```
Require R
Usage:
R binomialtest.R --no-save input.bed output.bed 0.2 < binomialtest.R
```
Description:

This function performs the binomial test to determine the signicant accumulation (20% of reads) of reads at certain 3'ends.





