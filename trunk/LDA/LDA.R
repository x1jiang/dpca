lap<-function(x, n)#n Laplacian noise whose scales are x
{
	z<-runif(n)
	y<-log(runif(n)) * ((z > 0.5) * 2 - 1)
	return (x * y)
}
main<-function()
{
	#load data
	x0<-read.table(inputfile)
	cat<-x0[1,]
	p<-dim(x0)[2]-1
	n<-dim(x0)[1]-1
	label<-(x0[2:(n+1),p+1]=="1")
	x<-x0[2:(n+1),1:p]
	nom<-(cat=='N')#whether nominal or not
	#count the total variables and find the names
	vc<-0
	vcl<-c()
	hash<-list()
	for (i in 1:p)
	{
		if (nom[i])
		{
			le<-unique(x[,i])
			vcl<-c(vcl,vc)
			vc<-vc+length(le)
			hash<-c(hash,list(le))
		}
		else
		{
			vcl<-c(vcl,vc)
			vc<-vc+1
			hash<-c(hash,list(""))
		}
	}
	vcl<-c(vcl,vc)
	vcl<-vcl+1
	data<-matrix(0,n,vc)# the data matrix
	minv<- rep(0,p)
	maxv<-rep(0,p)
	#step1, transform, normalize and split the data into two classes
	for (i in 1:p)
	{
		if (nom[i])
		{
			v <- x[,i]
			ov <- order(v)
			v<-v[ov]
			j<-1
			rec<-v[1]
			for (k in 1:n)
			{
				if (v[k]==rec)
				{
					data[ov[k],vcl[i]-1+j] <- 1
				}
				else
				{
					j<-j+1
					rec<-v[k]
					data[ov[k],vcl[i]-1+j] <- 1
				}
			}
		}
		else
		{
			v<-x[,i]
			v<-as.numeric(matrix(unlist(v),length(v),1))
			minv[i]<-min(v)
			maxv[i]<-max(v)
			data[,vcl[i]]<-(v-minv[i])/(maxv[i]-minv[i])
		}
	}
	data1<-data[label,]
	data2<-data[!label,]
	n1<- sum(label)
	n2<-n-n1
	#step2, generate noisy statistics
	p2<-sum(nom)-1
	p1<-p-p2
	s<-(2*p1+p2*2+((p+1)*p+p2*(p2+1))/2)/epsilon
	if (n1 != 0)
  {
  	c1m1<-rep(1/n1,n1)%*%data1
  	c1m1<-c1m1+lap(s/n1,vc)
  }
  if (n1 == 0)
  {
    c1m1<-rep(0,vc)
  }
  if (n2 != 0)
  {
  	c2m1<-rep(1/n2,n2)%*%data2
  	c2m1<-c2m1+lap(s/n2,vc)
  }
  if (n2 == 0)
  {
    c2m1<-rep(0,vc)
  }
	m2<-t(data)%*%data/n

	
	c2m1<-c2m1+lap(s/n2,vc)
	for (i in 1:vc)
	{
		m2[i,i:vc]<-m2[i,i:vc]+lap(s/n,vc-i+1)
	}	
	for (i in 1:vc)
	{
		m2[i,1:(i-1)]<-m2[1:(i-1),i]
	}
	
	#step3, sample new samples with C program since it is computation-intensive
	#reorder the parameters for convenience in sampling
	ords<-c()
	vars<-c()
	cons<-c()
	for (i in 1:p)
	{
		if (nom[i])
		{
			cons<-c(cons,vcl[i+1]-vcl[i])
			vars<-c(vars,-1)
		}
		else
		{
      j<-vcl[i]
			vars<-c(vars,max(m2[j,j],0))
			cons<-c(cons,-1)
		}
	}
	vars<-order(vars)[(p-p1+1):p]
  if (p2 != 0)
  {
  	cons2<-order(cons)[(p-p2+1):p]
  	for (i in 1:p2)
  	{
  	  ords<-c(ords,vcl[cons2[i]]:(vcl[cons2[i]+1]-1))
  	}
  	cons<-cons[cons2]
  }
  if (p2 == 0)
  {
    cons2<-c()
    cons<-c()
  }

	ords<-c(ords,vcl[vars])
	c1m1<-c1m1[ords]
	c2m1<-c2m1[ords]
	m2<-m2[ords,ords]
	#write to the file
	cons<-c(cons,rep(0,vc-length(cons)))
	outdata<-rbind(rep(vc,vc),rep(n1,vc), rep(n2,vc),cons, c1m1, c2m1, m2)
	write.table(outdata,"cppin.txt",sep=" ", col.names=FALSE,row.names=FALSE)
	#use cpp to sample
	system2("sample.bat","",TRUE)
	#read from file
	data<-read.table("cppout.txt")
	data<-matrix(unlist(data),dim(data)[1],dim(data)[2])
	#restore the order
	data<-data[,c(order(c(cons2,vars)),dim(data)[2])]

	#transform to the original form and save in y
	y<-x0[1:n,]
	for (i in 1:p)
	{
		if (nom[i])
		{
			labels<-unlist(hash[i])
			y[,i]<-labels[data[,i]-min(data[,i])+1]
		}
		else
		{
			y[,i]<-data[,i]*(maxv[i]-minv[i])+minv[i]
		}
	}
  y[,p+1]<-data[,p+1]
	#write the output
	cat<-matrix(unlist(cat),1,p+1)
	write.table(cat, outputfile, col.names=colnames(x0), row.names=FALSE, sep="")
  write.table(y, outputfile, append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ")
}
#main<-function(inputfile,outputfile,epsilon,k)
#{
#
#	data1<-cbind(LDAoneClass(x1,p,n1,nom),"1")
#	data2<-cbind(LDAoneClass(x2,p,n2,nom),"0")
#	write.table(rbind(cat,data1,data2), outputfile, col.names=colnames(x), row.names=FALSE, sep=" ")
#}

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
	if (c=="-h")
		hl=TRUE
}
if (epsilon <= 0)
	epsilon <- 1
if (!hl)
{
	main()
}
if (hl)
{
	help()
}
