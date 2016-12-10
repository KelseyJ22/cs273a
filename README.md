Epigenomic Data Acquisition
- - - - - - - - - - - - - - 
1. We downloaded h3k4m31.wig, h2k4me3.wig, h3k9ac.wig, h3k27me3.wig, and h3k36me3.wig from the Roadmap Epigenomics Project, as directed.
2. We wrote a script (fasta_to_bed.py) to convert the fasta-ike files storing Vista Enhancer Browser data to a .bed file we could use with the UCSC utilities.
3. The binaries wigToBigWig and bigWigAverageOverBed are downloaded from the UCSC Genome Browser and used as "wigToBigWig input.wig chrom.sizes myBigWig.bw" (with input.wig as zipped or unzipped wig file and chrom.sizes computed by the fetchChromSizes utility) and "./bigWigAverageOverBed in.bw in.bed out.tab", respectively. We also downloaded and used fetChromosomes to get the necessary chromosome information.
4. That resulted in the files h3k4m31.tab, h2k4me3.tab, h3k9ac.tab, h3k27me3.tab, and h3k36me3.tab with data of the form 

0	639	639	33	0.0516432	0.0516432

1	1641	1641	1116	0.680073	0.680073

2	1309	1309	317	0.24217	0.24217

3	1162	1162	200	0.172117	0.172117

4	1390	1390	357	0.256835	0.256835

5	1568	1568	450	0.28699	0.28699

6	1599	1599	440	0.275172	0.275172


etc

5. We then wrote a parser (map_avg_to_seq.py) to compare the .tab files to the .bed file and connect results generated by bigWigAverageOverBed to a sequence ftom the Vista Enhancer Browser. 
