# creation of default
sample_f<-read.table("sample_f.csv", header=T, sep = ",", row.names=1)
GSE58137.samples<-read.table("GSE58137.samples.txt", header=T, sep = "\t", row.names=1)
GSE58137.entrez <- data.frame(fread("GSE58137.entrez.pcl", header = T,sep="\t"), row.names = 1)
blood.f <- data.frame(fread("blood.f.pcl", header = T,sep="\t"), row.names = 1)
save.image("default.RData")