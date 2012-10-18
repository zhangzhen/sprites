x <- read.table("chr22_report.txt",sep="\t",header=TRUE)
x1 <- x[x$class=="DEL",]
nrow(x1)
y <- x1[,3]-x1[,2]
max(y)
length(y[y>=50])
hist(y,50)
