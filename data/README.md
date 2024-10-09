# how to download
Download the preprocessed data by copy-past the following to any web browser and download the zip file.

https://cbb.ittc.ku.edu/software/kmer/data.zip


# how to use
After downloading, you should unzip the data.zip and put all included files under the data/ folder. 
You should see 4 files included:
1. **otu.count.csv**
  - the **comma-separated values** file of the raw read count of the taxonomic profile. 
  - Rows are species-level taxa and columns are samples. 
1. **pfam.count.tsv** 
  - the **tab-separated-separated values** file of the raw read count of the Pfam profile. 
  - Rows are protein families and columns are samples. 
1. **k30.count.tsv**
  - the **tab-separated-separated values** file of the kmer count profile. 
  - Rows are kmer ids and columns are samples. 
  - you can use the tellkmer.py function to find out the actual kmer sequences.
1. **pfam.acclengname.csv**
  - the comma-separated values file of information regarding Pfam features, which is used to normalize the raw read count to RPKM. 


