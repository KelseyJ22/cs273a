Vista Enhancer Browser Data Acquisition
- - - - - - - - - - - - - -
Enhancer sequence data was downloaded from the VISTA Enhancer Browser at https://enhancer.lbl.gov/. The data can be downloaded under the 'Experimental Data' tab, with the 'Download Data' link at the top. The data is a list of sequences with labels for whether the sequence is an enhancer (positive or negative) as well as the regions in which the enhancer is active if it is postiive. The data is then parsed using project.py (run with command 'python project.py'). The parser assumes the VISTA data is named 'vistaEnhancerBrowser.txt'. It will save the features (including epigenomic data) in a .pkl file named 'pickled_vista_data_k_*.pkl'. The value of k (in the k-mers) can be edited in the code. The epigenomic data (.tab files) should be created before running the parser (see below for instructions). The .pkl file will contain all the feature vectors and enhancer type labels for each sequence in the VISTA data (in the pandas dataframe format). For k = 6, the parser should run in about 15 minutes. For k = 5, it should run in about 3 minutes. For k = 4, it should run in around 45 seconds. For k = 3, it should run in around 15 seconds.


Using Beer Lab's Data
- - - - - - - - - - - - - -
To train and run the classifiers on data from the Beer lab (the lab in the original paper), you should download the data from http://www.beerlab.org/p300enhancer/ (under the Training Data Sets). These files are labeled with the experiments that were run to obtain the data, which are described in the paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3227105/). They are FASTA format files of sequences that are shown to be enhancers for the appropriate region. Some of the files are labeled 'random*.fa'. These files are negative data - they are randomly generated sequences that are not enhancers. The number denotes the number of samples in the negative data set. When training classifiers, an approximately equal number of positive and negative data samples should be used (the random files can be shortened if necessary). In addition, you should select the random data set which has approximately equal distribution of sequence lengths as the positive data set (this can be done just by observing the .fa files).

To parse the data, run beer_data_parser.py (with command 'python beer_data_parser.py'). This parser does the same thing as the VISTA data parser. However, the positive data file and negative data file must both be specified (as pos_file and neg_file in the main() function). In addition, the output .pkl file should be edited in to the out_file variable; the parser will not automatically generate the correct filename.

The parser runs in reasonable time on the smaller data sets (< 5000 total positive and negative samples). However, the larger data sets will likely take an unreasonable amount of time to parse.


Training and Running the Classifiers
- - - - - - - - - - - - - -
The classifiers can be trained and run using the plot_*_classifier.py scripts (run with command 'python plot_*_classifier.py'). The classifiers are KNN, naive Bayes, logistic regression, multilayer perceptron, random forest, and SVM. They will perform binary classification for whether a ssequence is an enhancer active in a specified region (or generally). For these scripts, the 'filename' variable in the 'read_pickle()' function should be changed to the .pkl data file for the parsed data that you want to use for the classifier. In addition, the code needs to be edited to correspond to the type of data used to generate the feature vectors and labels. If VISTA data was used, LABEL_SIZE should be set to 8. enhancer_labels should be set to labels_matrix[labels_matrix.columns[x]], where x is a number from 0 to 7 corresponding to the label you want the classifier to classify on (the order of the labels is given in the comments just above the variable declaration). If Beer lab data was used, LABEL_SIZE should be set to 1. enhancer_labels must be set to labels_matrix[labels_matrix.columns[0]], since the label is specified by the data set used (e.g. enh_fb.fa for forebrain enhancer classification). 

The scripts will print out statistics on the classifier performance, as well as saving an image of the ROC curve for the classifier. Depending on the data used and the value of k used to generate the features, the classifier may take a couple of minutes to run.


Epigenomic Data Acquisition
- - - - - - - - - - - - - - 
1. We downloaded h3k4me1.wig, h2k4me3.wig, h3k9ac.wig, h3k27me3.wig, and h3k36me3.wig from the Roadmap Epigenomics Project, as directed.
2. We wrote a script (fasta_to_bed.py) to convert the fasta-ike files storing Vista Enhancer Browser data to a .bed file we could use with the UCSC utilities.
3. The binaries wigToBigWig and bigWigAverageOverBed are downloaded from the UCSC Genome Browser and used as "wigToBigWig input.wig chrom.sizes myBigWig.bw" (with input.wig as zipped or unzipped wig file and chrom.sizes computed by the fetchChromSizes utility) and "./bigWigAverageOverBed in.bw in.bed out.tab", respectively. We also downloaded and used fetChromosomes to get the necessary chromosome information.
4. That resulted in the files H3K4me1.tab, H2K4me3.tab, H3K9ac.tab, H3K27me3.tab, and H3K36me3.tab with data of the form 

0	639	639	33	0.0516432	0.0516432

1	1641	1641	1116	0.680073	0.680073

2	1309	1309	317	0.24217	0.24217

3	1162	1162	200	0.172117	0.172117

4	1390	1390	357	0.256835	0.256835

5	1568	1568	450	0.28699	0.28699

6	1599	1599	440	0.275172	0.275172


etc

5. We then wrote a parser (map_avg_to_seq.py) to compare the .tab files to the .bed file and connect results generated by bigWigAverageOverBed to a sequence from the Vista Enhancer Browser. 
