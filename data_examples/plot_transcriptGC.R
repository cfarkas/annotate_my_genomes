library(car)
data1<-read.table(file="GC_nuc", header=F, sep = "\t")
data2 <- data1[order(data1$V12),]
head(data2)
pdf("transcript_GC_content.pdf", height=11, width=11)
scatterplot(data2$V12 ~ data2$V5, data = data2, xlab="% GC content", ylab= "Transcript Length (bp)", col="darkgreen")
dev.off()