lap<-function(x, n)#n Laplacian noise whose scales are x
{
	z<-runif(n)
	y<-log(runif(n)) * ((z > 0.5) * 2 - 1)
	return (x * y)
}
main<-function(inputfile,outputfile,epsilon,k)
{
	#load data
	x<-read.table(inputfile)
	cat<-x[1,]
	p<-dim(x)[2]
	n<-dim(x)[1]-1
	nom<-(cat=='N')#whether nominal or not
	
	#count the total variables and find the names
	vc<-0
	vc1<-c()
	hash<-list()
	for (i in 1:p)
	{
		if (nom[i])
		{
			le<-unique(x[2:(n+1),i])
			vc1<-c(vc1,vc)
			vc<-vc+length(le)
			hash<-c(hash,list(le))
		}
		else
		{
			vc1<-c(vc1,vc)
			vc<-vc+1
			hash<-c(hash,list(""))
		}
	}
	vc1<-c(vc1,vc)
	vc1<-vc1+1
	data<-matrix(0,n,vc)# the data matrix
	minv<- rep(0,p)
	maxv<-rep(0,p)
	#step1, transform, normalize and store the data
	for (i in 1:p)
	{
		if (nom[i])
		{
			v <- x[2:(n+1),i]
			ov <- order(v)
			v<-v[ov]
			j<-1
			rec<-v[1]
			for (k in 1:n)
			{
				if (v[k]==rec)
				{
					data[ov[k],vc1[i]-1+j] <- 1
				}
				else
				{
					j<-j+1
					rec<-v[k]
					data[ov[k],vc1[i]-1+j] <- 1
				}
			}
		}
		else
		{
			v<-x[2:(n+1),i]
			v<-as.numeric(matrix(unlist(v),length(v),1))
			minv[i]<-min(v)
			maxv[i]<-max(v)
			data[,vc1[i]]<-(v-minv[i])/(maxv[i]-minv[i])
		}
	}

	#step2, get the statistis
	p2<-sum(nom)
	p1<-p-p2
	m1<-rep(1/n,n)%*%data
	s<-2*(p1+p2*2+((p+1)*p+(p2+1)*p2)/2)/epsilon
	m2<-t(data)%*%data/n

	#step3, generate noisy statistics
	m1<-m1+lap(s/n,vc)
	for (i in 1:vc)
	{
		m2[i,i:vc]<-m2[i,i:vc]+lap(s/n,vc-i+1)
	}	
	for (i in 1:vc)
	{
		m2[i,1:(i-1)]<-m2[1:(i-1),i]
	}
	#step4, get uk
	covar<-m2-t(m1)%*%m1

	eigens<-eigen(covar)
	if (k<1)
		k<-1
	if (k>vc)
		k<-vc
	u<-matrix(unlist(eigens$vectors),vc,vc)
	uk<-u[,1:k]
	
	#step5, get noisy data
	noisescale<-2*sqrt(k*vc)/epsilon
	md<-(data-rep(1/n,n)%*%m1)%*%uk+lap(noisescale, n*vc)
	md<-md%*%t(uk)
	for (i in 1:vc)
		data[,i]<-md[,i]+m1[,i]
	
	#change form, select one value for each nominal variable
	for (i in 1:p)
	{
		if (nom[i])
		{
			range<-vc1[i]:(vc1[i+1]-1)
			for (j in 1:n)
			{
				l<-data[j,range]
				l<-(max(l)==l)
				data[j,range]<-l
			}
		}
	}

	#transform to the original form and save in y
	y<-x[1:n,1:p]
	for (i in 1:p)
	{
		if (nom[i])
		{
			labels<-unlist(hash[i])
			for (j in vc1[i]:(vc1[i+1]-1))
			{
				label<-(data[,j]==1)
				y[label,i]<-labels[j-vc1[i]+1]
			}
		}
		else
		{
			y[,i]<-data[,vc1[i]]*(maxv[i]-minv[i])+minv[i]
		}
	}
	#write the output
	cat<-matrix(unlist(cat),1,p)
  write.table(cat, outputfile, col.names=colnames(x), row.names=FALSE, sep="")
	write.table(y, outputfile, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ")
}
help<-function()
{
	print("Usage: PCA.bat [--options] [args]")
	print(" ")
	print("--options accepted are")
	print("  -h              Print help information and exit")
	print(" ")
	print("args accepted are")
	print("  -e              Value of Privacy budget, epsilon. epsilon=1 by default")
	print("  -k              The number of eigen vectors used in PCA. It will be changed to 1 if less than 1, and the number of variables if it is larger than that. k=1 by default")
	print("  -o              The file to which the noisy data are written")
	print("  -i              The file from which the original data are read")
	print(" ")
	print("Data format:")
	print("Each data set shall be a matrix. The first line is the name of variables, while the second line indicates whether this is nominal or continuous using labels N and C. Nominal variables will be changed to dummy variables. The output file will be in the same format.")
}


#read inputfile, outputfile and epsilon from arguments
argv<-commandArgs(TRUE)
argc<-length(argv)
i<-1
k<-1#number of components used
outputfile<-"output.txt"
inputfile<-"input.txt"
epsilon<-1
hl<-FALSE
while (i <= argc)
{
	c <- argv[i]
	i<-i+1
	if (c =="-o")
		outputfile<-argv[i]
	if (c=="-i")
		inputfile<-argv[i]
	if (c=="-e")
		epsilon<-as.numeric(argv[i])
	if (c=="-k")
		k<-as.numeric(argv[i])
	if (c=="-h")
		hl=TRUE
}
if (epsilon <= 0)
	epsilon <- 1
if (!hl)
{
	main(inputfile,outputfile,epsilon,k)
}
if (hl)
{
	help()
}
