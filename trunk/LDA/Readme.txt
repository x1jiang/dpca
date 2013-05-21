This is a differentially private LDA program in R.

The data set format is as follows: it shall be a matrix whose first line is the name of variables while the second line indicates whether this variable is nominal (N) or continuous (C). The output data set has the same format.

Since LDA is a classification algorithm, the nominal response variable is in the last column. The response variable can be either 1 or 0.

There are one option, "-h" to output the help information.

There are 4 parameters that can be set: names of input and output files, the privacy budget epsilon. A sample is:
/Path/to/Rscript.exe LDA.R -i samplein.txt -o sampleout.txt -e 1