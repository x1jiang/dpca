This is a differentially private PCA program in R.

The data set format is as follows: it shall be a matrix whose first line is the name of variables while the second line indicates whether this variable is nominal (N) or continuous (C). The output data set has the same format.

There are one option, "-h" to output the help information.

There are 4 parameters that can be set: names of input and output files, the privacy budget epsilon and the number of components used, k. A sample is:
/Path/to/Rscript.exe PCA.R -i samplein.txt -o sampleout.txt -e 1 -k 10