# VicRcodes-mutation_bars


This code was designed by Victor Jimenez Vasquez - vr.jimenez.vs@gmail.com.

![graphic1](https://user-images.githubusercontent.com/89874227/156581714-857a724b-d9fd-4abf-a4df-84b5e8c1cfda.jpg)

## Intro

Genomic surveillance of SARS-CoV-2 include the identification of variable sites at aminoacid and genome levels. Mutation_bars allows the user the identification of variable sites in a given data set composed of FASTA files, metadata downloaded from GISAID and lineage report obtained from PANGOLIN. 

## Usage 
```r
#1# setwd and load libraries#
setwd("F:/sarscov2/variations/phylo/progress_2020/origins/filogenia_2020")

data.frame(dir())
library(seqinr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("mutation_bars.R")

#2# identify the first element containing the first coding (Envelope) fasta file (x) and the last element containing the last coding (Spike) fasta file (y)#
data.frame(dir())
r <- data.frame(dir()[x:y])

#3# read fasta files#
E <- read.fasta(r[1,])
M <- read.fasta(r[2,])
N <- read.fasta(r[3,])
orf1a <- read.fasta(r[4,])
orf1b <- read.fasta(r[5,])
orf3a <- read.fasta(r[6,])
orf6 <- read.fasta(r[7,])
orf7a <- read.fasta(r[8,])
orf7b <- read.fasta(r[9,])
orf8 <- read.fasta(r[10,])
orf9b <- read.fasta(r[11,])
s <- read.fasta(r[12,])

#4# read pango file#
pango <- read.csv("lineage_report.csv", header=TRUE)
names(pango)
pango <- pango[,c(1:2)]
names(pango) <- c("taxon","lineage")

#5# read metadata#
meta <- read.csv("2020_metadata.tsv", header=TRUE, sep="\t")
names(meta)

#6# optional: you may select some columns#
meta <- meta[,c(1:7,9)]
names(meta) <- c("taxon","virus","gisaid_epi_isl","genbank_accession","date","region","country","location")

#7# transform each fasta to data.frame#
a <- fastameta_to_df(fasta=E,gene="E",pango=pango,label=2)
b <- fasta_to_df(fasta=M,gene="M")
c <- fasta_to_df(fasta=N,gene="N")
d <- fasta_to_df(fasta=orf1a,gene="ORF1a")
e <- fasta_to_df(fasta=orf1b,gene="ORF1b")
f <- fasta_to_df(fasta=orf3a,gene="ORF3a")
g <- fasta_to_df(fasta=orf6,gene="ORF6")
h <- fasta_to_df(fasta=orf7a,gene="ORF7a")
i <- fasta_to_df(fasta=orf7b,gene="ORF7b")
j <- fasta_to_df(fasta=orf8,gene="ORF8")
k <- fasta_to_df(fasta=orf9b,gene="ORF9b")
l <- fasta_to_df(fasta=s,gene="S")

#8# merge data.frames by taxon column#
m <- merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(a,b,by="taxon", all.x=TRUE),
c,by="taxon", all.x=TRUE),
d,by="taxon", all.x=TRUE),
e,by="taxon", all.x=TRUE),
f,by="taxon", all.x=TRUE),
g,by="taxon", all.x=TRUE),
h,by="taxon", all.x=TRUE),
i,by="taxon", all.x=TRUE),
j,by="taxon", all.x=TRUE),
k,by="taxon", all.x=TRUE),
l,by="taxon", all.x=TRUE)

#9# estimate the number of rows and columns to your data#
dim(m)
names(m)[1:20]

#10# generate final input#
n <- merge(meta,m,by="taxon", all.x=FALSE)
dim(n)
names(n)[1:100]

#11# final command#
res <- barras2(data=n,linaje="BA.1",genomas=10,run="peru_2020",label="coding",inic=10)

#12# arguments#
data : input containing sites for all coding regions and metadata
linaje : lineage
genomas : the minimum number of genomes in which the program will search for variable sites in your input
run : word to be print in the resulting barplot
label : output file name 
inic : the first column number containing aminoacid information in your input

```

## Output
1. A barplot for the aminoacid frequency per variable site identified in a minimum number of genomes as especified in "genomas" argument. 
2. A data table in csv format containing a metadata plus variable sites columns identified in a minimum number of genomes as especified in "genomas" argument. 
 
## Usage
1. Use "mutation_bars" to obtain a graphical representation of the variable sites represented by a minimum number of genomes in your data
2. Use "mutation_bars" to obtain a complete metadata in csv format including mutation sites that you can map in your phylogenetic tree for a visual inspection in https://microreact.org/ 

This code was designed by Victor Jimenez Vasquez - vr.jimenez.vs@gmail.com 

Vic
