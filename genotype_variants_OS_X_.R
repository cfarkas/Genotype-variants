#!/bin/bash

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(tcltk2)

WT_variants <- file.choose()
WT<-read.table(WT_variants, skip="##")
WT<-subset(WT, select = c(1, 2))
head(WT)
dim(WT)
Null_variants <- file.choose()
Null<-read.table(Null_variants, skip="##")
Null<-subset(Null, select = c(1, 2))
head(Null)
dim(Null)
########## subsetting chromosomes ###############

WT_chr1<-subset(WT, V1 == "chr1")
WT_chr2<-subset(WT, V1 == "chr2")
WT_chr3<-subset(WT, V1 == "chr3")
WT_chr4<-subset(WT, V1 == "chr4")
WT_chr5<-subset(WT, V1 == "chr5")
WT_chr6<-subset(WT, V1 == "chr6")
WT_chr7<-subset(WT, V1 == "chr7")
WT_chr8<-subset(WT, V1 == "chr8")
WT_chr9<-subset(WT, V1 == "chr9")
WT_chr10<-subset(WT, V1 == "chr10")
WT_chr11<-subset(WT, V1 == "chr11")
WT_chr12<-subset(WT, V1 == "chr12")
WT_chr13<-subset(WT, V1 == "chr13")
WT_chr14<-subset(WT, V1 == "chr14")
WT_chr15<-subset(WT, V1 == "chr15")
WT_chr16<-subset(WT, V1 == "chr16")
WT_chr17<-subset(WT, V1 == "chr17")
WT_chr18<-subset(WT, V1 == "chr18")
WT_chr19<-subset(WT, V1 == "chr19")
WT_chrX<-subset(WT, V1 == "chrX")

Null_chr1<-subset(Null, V1 == "chr1")
Null_chr2<-subset(Null, V1 == "chr2")
Null_chr3<-subset(Null, V1 == "chr3")
Null_chr4<-subset(Null, V1 == "chr4")
Null_chr5<-subset(Null, V1 == "chr5")
Null_chr6<-subset(Null, V1 == "chr6")
Null_chr7<-subset(Null, V1 == "chr7")
Null_chr8<-subset(Null, V1 == "chr8")
Null_chr9<-subset(Null, V1 == "chr9")
Null_chr10<-subset(Null, V1 == "chr10")
Null_chr11<-subset(Null, V1 == "chr11")
Null_chr12<-subset(Null, V1 == "chr12")
Null_chr13<-subset(Null, V1 == "chr13")
Null_chr14<-subset(Null, V1 == "chr14")
Null_chr15<-subset(Null, V1 == "chr15")
Null_chr16<-subset(Null, V1 == "chr16")
Null_chr17<-subset(Null, V1 == "chr17")
Null_chr18<-subset(Null, V1 == "chr18")
Null_chr19<-subset(Null, V1 == "chr19")
Null_chrX<-subset(Null, V1 == "chrX")

########## subsetting coordinates ###############

WT_chr1<-subset(WT_chr1, select = c(2))
WT_chr2<-subset(WT_chr2, select = c(2))
WT_chr3<-subset(WT_chr3, select = c(2))
WT_chr4<-subset(WT_chr4, select = c(2))
WT_chr5<-subset(WT_chr5, select = c(2))
WT_chr6<-subset(WT_chr6, select = c(2))
WT_chr7<-subset(WT_chr7, select = c(2))
WT_chr8<-subset(WT_chr8, select = c(2))
WT_chr9<-subset(WT_chr9, select = c(2))
WT_chr10<-subset(WT_chr10, select = c(2))
WT_chr11<-subset(WT_chr11, select = c(2))
WT_chr12<-subset(WT_chr12, select = c(2))
WT_chr13<-subset(WT_chr13, select = c(2))
WT_chr14<-subset(WT_chr14, select = c(2))
WT_chr15<-subset(WT_chr15, select = c(2))
WT_chr16<-subset(WT_chr16, select = c(2))
WT_chr17<-subset(WT_chr17, select = c(2))
WT_chr18<-subset(WT_chr18, select = c(2))
WT_chr19<-subset(WT_chr19, select = c(2))
WT_chrX<-subset(WT_chrX, select = c(2))

Null_chr1<-subset(Null_chr1, select = c(2))
Null_chr2<-subset(Null_chr2, select = c(2))
Null_chr3<-subset(Null_chr3, select = c(2))
Null_chr4<-subset(Null_chr4, select = c(2))
Null_chr5<-subset(Null_chr5, select = c(2))
Null_chr6<-subset(Null_chr6, select = c(2))
Null_chr7<-subset(Null_chr7, select = c(2))
Null_chr8<-subset(Null_chr8, select = c(2))
Null_chr9<-subset(Null_chr9, select = c(2))
Null_chr10<-subset(Null_chr10, select = c(2))
Null_chr11<-subset(Null_chr11, select = c(2))
Null_chr12<-subset(Null_chr12, select = c(2))
Null_chr13<-subset(Null_chr13, select = c(2))
Null_chr14<-subset(Null_chr14, select = c(2))
Null_chr15<-subset(Null_chr15, select = c(2))
Null_chr16<-subset(Null_chr16, select = c(2))
Null_chr17<-subset(Null_chr17, select = c(2))
Null_chr18<-subset(Null_chr18, select = c(2))
Null_chr19<-subset(Null_chr19, select = c(2))
Null_chrX<-subset(Null_chrX, select = c(2))

par(mfrow=c(2,10))

WT_chr1$Genotype<-"WT"
Null_chr1$Genotype<-"Null"
Chr1 <- rbind(WT_chr1, Null_chr1)
p1 <- ggplot(Chr1, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr1") + ylab("Variants") +labs(title="Chr1") + theme(legend.position="none")

WT_chr2$Genotype <- "WT"
Null_chr2$Genotype <- "Null"
Chr2 <- rbind(WT_chr2, Null_chr2)
p2 <- ggplot(Chr2, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr2") + ylab("Variants") +labs(title="Chr2") + theme(legend.position="none")

WT_chr3$Genotype <- "WT"
Null_chr3$Genotype <- "Null"
Chr3 <- rbind(WT_chr3, Null_chr3)
p3 <- ggplot(Chr3, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr3") + ylab("Variants") +labs(title="Chr3") + theme(legend.position="none")

WT_chr4$Genotype<-"WT"
Null_chr4$Genotype<-"Null"
Chr4 <- rbind(WT_chr4, Null_chr4)
p4 <- ggplot(Chr4, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr4") + ylab("Variants") +labs(title="Chr4") + theme(legend.position="none")

WT_chr5$Genotype<-"WT"
Null_chr5$Genotype<-"Null"
Chr5 <- rbind(WT_chr5, Null_chr5)
p5 <- ggplot(Chr5, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr5") + ylab("Variants") +labs(title="Chr5") + theme(legend.position="none")

WT_chr6$Genotype<-"WT"
Null_chr6$Genotype<-"Null"
Chr6 <- rbind(WT_chr6, Null_chr6)
p6 <- ggplot(Chr6, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr6") + ylab("Variants") +labs(title="Chr6") + theme(legend.position="none")

WT_chr7$Genotype<-"WT"
Null_chr7$Genotype<-"Null"
Chr7 <- rbind(WT_chr7, Null_chr7)
p7 <- ggplot(Chr7, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr7") + ylab("Variants") +labs(title="Chr7") + theme(legend.position="none")

WT_chr8$Genotype<-"WT"
Null_chr8$Genotype<-"Null"
Chr8 <- rbind(WT_chr8, Null_chr8)
p8 <- ggplot(Chr8, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr8") + ylab("Variants") +labs(title="Chr8") + theme(legend.position="none")

WT_chr9$Genotype<-"WT"
Null_chr9$Genotype<-"Null"
Chr9 <- rbind(WT_chr9, Null_chr9)
p9 <- ggplot(Chr9, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr9") + ylab("Variants") +labs(title="Chr9") + theme(legend.position="none")

WT_chr10$Genotype<-"WT"
Null_chr10$Genotype<-"Null"
Chr10 <- rbind(WT_chr10, Null_chr10)
p10 <- ggplot(Chr10, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr10") + ylab("Variants") +labs(title="Chr10") 

WT_chr11$Genotype<-"WT"
Null_chr11$Genotype<-"Null"
Chr11 <- rbind(WT_chr11, Null_chr11)
p11 <- ggplot(Chr11, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr11") + ylab("Variants") +labs(title="Chr11") + theme(legend.position="none")

WT_chr12$Genotype<-"WT"
Null_chr12$Genotype<-"Null"
Chr12 <- rbind(WT_chr12, Null_chr12)
p12 <- ggplot(Chr12, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr12") + ylab("Variants") +labs(title="Chr12") + theme(legend.position="none")

WT_chr13$Genotype<-"WT"
Null_chr13$Genotype<-"Null"
Chr13 <- rbind(WT_chr13, Null_chr13)
p13 <- ggplot(Chr13, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr13") + ylab("Variants") +labs(title="Chr13") + theme(legend.position="none")

WT_chr14$Genotype<-"WT"
Null_chr14$Genotype<-"Null"
Chr14 <- rbind(WT_chr14, Null_chr14)
p14 <- ggplot(Chr14, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr14") + ylab("Variants") +labs(title="Chr14") + theme(legend.position="none")

WT_chr15$Genotype<-"WT"
Null_chr15$Genotype<-"Null"
Chr15 <- rbind(WT_chr15, Null_chr15)
p15 <- ggplot(Chr15, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr15") + ylab("Variants") +labs(title="Chr15") + theme(legend.position="none")

WT_chr16$Genotype<-"WT"
Null_chr16$Genotype<-"Null"
Chr16 <- rbind(WT_chr16, Null_chr16)
p16 <- ggplot(Chr16, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr16") + ylab("Variants") +labs(title="Chr16") + theme(legend.position="none")

WT_chr17$Genotype<-"WT"
Null_chr17$Genotype<-"Null"
Chr17 <- rbind(WT_chr17, Null_chr17)
p17 <- ggplot(Chr17, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr17") + ylab("Variants") +labs(title="Chr17") + theme(legend.position="none")

WT_chr18$Genotype<-"WT"
Null_chr18$Genotype<-"Null"
Chr18 <- rbind(WT_chr18, Null_chr18)
p18 <- ggplot(Chr18, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr18") + ylab("Variants") +labs(title="Chr18") + theme(legend.position="none")

WT_chr19$Genotype<-"WT"
Null_chr19$Genotype<-"Null"
Chr19 <- rbind(WT_chr19, Null_chr19)
p19 <- ggplot(Chr19, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr19") + ylab("Variants") +labs(title="Chr19") + theme(legend.position="none")

WT_chrX$Genotype<-"WT"
Null_chrX$Genotype<-"Null"
ChrX <- rbind(WT_chrX, Null_chrX)
pX <- ggplot(ChrX, aes(V2, fill = Genotype)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrX") + ylab("Variants") +labs(title="ChrX")



grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}

dsamp <- diamonds[sample(nrow(diamonds), 1000), ]
g <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 2, nrow=5)
ggsave("graph1.pdf", g, width = 25, height = 50, units = "cm")
g <- grid_arrange_shared_legend(p11, p12, p13, p14, p15, p16, p17, p18, p19, pX, ncol = 2, nrow=5)
ggsave("graph2.pdf", g, width = 25, height = 50, units = "cm")

### statistical tests ###

par(mar = rep(2, 4))
WTchr1<-hist(WT_chr1$V2)
Nullchr1<-hist(Null_chr1$V2)
a<-var.test(WTchr1$counts, Nullchr1$counts)
Chr1=a

WTchr2<-hist(WT_chr2$V2)
Nullchr2<-hist(Null_chr2$V2)
b<-var.test(WTchr2$counts, Nullchr2$counts)
Chr2=b

WTchr3<-hist(WT_chr3$V2)
Nullchr3<-hist(Null_chr3$V2)
c<-var.test(WTchr3$counts, Nullchr3$counts)
Chr3=c

WTchr4<-hist(WT_chr4$V2)
Nullchr4<-hist(Null_chr4$V2)
d<-var.test(WTchr4$counts, Nullchr4$counts)
Chr4=d

WTchr5<-hist(WT_chr5$V2)
Nullchr5<-hist(Null_chr5$V2)
e<-var.test(WTchr5$counts, Nullchr5$counts)
Chr5=e

WTchr6<-hist(WT_chr6$V2)
Nullchr6<-hist(Null_chr6$V2)
f<-var.test(WTchr6$counts, Nullchr6$counts)
Chr6=f

WTchr7<-hist(WT_chr7$V2)
Nullchr7<-hist(Null_chr7$V2)
g<-var.test(WTchr7$counts, Nullchr7$counts)
Chr7=g

WTchr8<-hist(WT_chr8$V2)
Nullchr8<-hist(Null_chr8$V2)
h<-var.test(WTchr8$counts, Nullchr8$counts)
Chr8=h

WTchr9<-hist(WT_chr9$V2)
Nullchr9<-hist(Null_chr9$V2)
i<-var.test(WTchr9$counts, Nullchr9$counts)
Chr9=i

WTchr10<-hist(WT_chr10$V2)
Nullchr10<-hist(Null_chr10$V2)
j<-var.test(WTchr10$counts, Nullchr10$counts)
Chr10=j

WTchr11<-hist(WT_chr11$V2)
Nullchr11<-hist(Null_chr11$V2)
k<-var.test(WTchr11$counts, Nullchr11$counts)
Chr11=k

WTchr12<-hist(WT_chr12$V2)
Nullchr12<-hist(Null_chr12$V2)
l<-var.test(WTchr12$counts, Nullchr12$counts)
Chr12=l

WTchr13<-hist(WT_chr13$V2)
Nullchr13<-hist(Null_chr13$V2)
m<-var.test(WTchr13$counts, Nullchr13$counts)
Chr12=m

WTchr14<-hist(WT_chr14$V2)
Nullchr14<-hist(Null_chr14$V2)
n<-var.test(WTchr14$counts, Nullchr14$counts)
Chr14=n

WTchr15<-hist(WT_chr15$V2)
Nullchr15<-hist(Null_chr15$V2)
o<-var.test(WTchr15$counts, Nullchr15$counts)
Chr15=o

WTchr16<-hist(WT_chr16$V2)
Nullchr16<-hist(Null_chr16$V2)
p<-var.test(WTchr16$counts, Nullchr16$counts)
Chr16=p

WTchr17<-hist(WT_chr17$V2)
Nullchr17<-hist(Null_chr17$V2)
q<-var.test(WTchr17$counts, Nullchr17$counts)
Chr17=q

WTchr18<-hist(WT_chr18$V2)
Nullchr18<-hist(Null_chr18$V2)
r<-var.test(WTchr18$counts, Nullchr18$counts)
Chr18=r

WTchr19<-hist(WT_chr19$V2)
Nullchr19<-hist(Null_chr19$V2)
s<-var.test(WTchr19$counts, Nullchr19$counts)
Chr19=s

WTchrX<-hist(WT_chrX$V2)
NullchrX<-hist(Null_chrX$V2)
t<-var.test(WTchrX$counts, NullchrX$counts)
Chr19=t

p.value<-c(a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value,g$p.value,h$p.value,i$p.value,j$p.value,k$p.value,l$p.value,m$p.value,n$p.value,o$p.value,p$p.value,q$p.value,r$p.value,s$p.value,t$p.value)
Chromosome<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX")
significance <- p.value < 0.05
matrix<-cbind(Chromosome, p.value, significance)
Summary<-data.frame(matrix)
Summary
pdf("Summary.pdf", height=11, width=8.5)
grid.table(Summary)
dev.off()
proc.time()
sessionInfo()