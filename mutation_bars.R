#########
### 1 ###
#########
#funcion para transformar una secuencia FASTA en un DATA.FRAME#
#este objeto data.frame puede guardarse como CSV#
#requiere el paquete "seqinr"#
#install.packages("seqinr")#
#el objeto "fasta" es un archivo fasta leido con la funcion read.fasta#

fasta_to_df <- function(fasta,gene){
gene = gene #el nombre de la region genica(E)#
b = fasta #archivo fasta leido con reas.fasta previamente#
n <- length(unlist(b[1], use.names=FALSE))
data <- data.frame(NA_col = rep(NA, n))
for(i in 1:length(b)) {
new_col <- unlist(b[i], use.names=FALSE)
data[ , i] <- new_col
colnames(data)[i] <- i
}
data2 <- t(data)
data3 <- data.frame(strain=names(b),data2)
names(data3) <- c("taxon",(paste0(rep(gene,(ncol(data3)-1)),"_",1:(ncol(data3)-1))))
data4 <- data3
}

#usage#
#fasta_to_df(fasta=fasta,gene=c("E","M","S","ORF1a")#

#########
### 2 ###
#########
#funcion para combinar una secuencia FASTA y output de PANGO en un DATA.FRAME#
#este objeto data.frame puede guardarse como CSV#
#requiere el paquete "seqinr"#

fastameta_to_df <- function(fasta,gene,pango,label){
label=label #etiqueta a incluir en el archivo salid de frecuencia linajes en "linajes_frecuencias_LABEL.csv"#
pango = pango #archivo resultante de la identificacion de linajes por el algoritmo pango#
gene = gene #el nombre de la region genica("E","M","S","ORF1a")#
b = fasta #archivo fasta leido con read.fasta previamente#

n <- length(unlist(b[1], use.names=FALSE))
data <- data.frame(NA_col = rep(NA, n))
for(i in 1:length(b)) {
new_col <- unlist(b[i], use.names=FALSE)
data[ , i] <- new_col
colnames(data)[i] <- i
}
data2 <- t(data)
data3 <- data.frame(strain=names(b),data2)
names(data3) <- c("taxon",(paste0(rep(gene,(ncol(data3)-1)),"_",1:(ncol(data3)-1))))
data4 <- data3

a <- merge(pango,data4, by="taxon", all=FALSE)
b <- data.frame(a, n= rep(1,dim(a)[1]))
lin <- aggregate(b$n, by=list(b$lineage), FUN=sum)
names(lin) <- c("linaje","freq")
ca <- lin[order(-lin$freq),]
lin2 <- lin[order(lin$linaje,lin$freq),]
lin2
title= paste0("linajes_frecuencias","_",label,".csv")
write.csv(ca, file=title, row.names=FALSE)
lin3 <- lin2
print(lin3)
print(a)
}

#usage#
#fastameta_to_df(fasta,region,pango,label)#

#########
### 3 ###
#########

#funcion para obtener plots de sitios variables por linaje#
#dichos sitios estan definidos si la variacion interna esta#
#presente en un numero minimo de genomas#

######## the function #######
#data = tab#
#linaje = "P.1.12"#
#genomas = 1#
#run = 31#
barras <- function(data,linaje,genomas,run,label){
data = data
linaje = linaje
run = run
genomas = genomas
names = unique(data$lineage)
ch <- data[data$lineage == linaje,]
ch2 <- ch[,14:(ncol(data)-1)]

#### WE WILL ESTIMATE THE VARIABLE SITES ####

ngenomes <- genomas
#prepare the loop#
freq <- 0
r <- matrix(nrow=dim(ch2)[2],ncol=1)
table <- matrix(nrow=dim(ch2)[2],ncol=1)
for(i in 1:dim(ch2)[2]) {
freq[i] <- as.data.frame(summary(as.factor(ch2[,i])))
table[i] <- data.frame(n=rep(i,length(unlist(freq[[i]]))),changes=unlist(freq[[i]]))
r[i] <- as.data.frame(names(summary(as.factor(ch2[,i]))))
print(freq[[i]])
print(r[[i]])
}
freq
table
r

##### WE WILL MERGE THE PRODUCED DATA #####
s <- data.frame(position=unlist(table),aa=unlist(r),freq=unlist(freq))
s
names(s)
s2 <- s[!s$aa=="x" & s$freq >= ngenomes,]
s2[order(-s2$freq),]
names(s2)
#write.csv(s, file="frequencies_S.csv",row.names=FALSE)#

##### WE WILL ESTIMATE THE FREQUENCIES FOR VARIABLE SITES WITH MORE THAN 2 STATES #####
s3 <- data.frame(s2, n=rep(1,dim(s2)[1]))
s3[1:25,]
freq <- aggregate(s3$n, by=list(s3$position), FUN=sum)
names(freq) <- c("position","freq")
freq2 <- freq[order(-freq$freq),]
names(freq2)
lists <- freq[freq$freq >= 2,]
lists[,1]
dim(lists)
length(lists[,1])
#write.csv(lists, file=paste0("sites_in",genomas,"genomes",".csv"),row.names=FALSE)#

S <- ch[order(ch$taxon),c(1:13,lists[,1]+13)]
names(S)
dim(S)
mutates_sites <- names(S)
write.csv(S, file=paste0(label,"_selected_positions.csv"), row.names=FALSE)

len <- length(lists[,1])

S2 <- data.frame(S[,c(14:dim(S)[2])],n=rep(1,dim(S)[1]))

names(S2)
print(S2)
v <- names(S[14])
nam <- function(len){
if(len == 1) return(c(v,"n"))
else if(len >= 2) return(names(S2))
}

b <- nam(len)
names(S2) <- b
S2

site <- c()
aa <- c()
freq <- c()

for (i in 1:(ncol(S2)-1)){
	AA <- aggregate(S2$n, by=list(S2[,i]), FUN=sum)[,1]
	BB <- aggregate(S2$n, by=list(S2[,i]), FUN=sum)[,2]
	site <- append(site,rep(colnames(S2)[i],length(AA)))
for (j in AA){
	aa <- append(aa, j)
}
for(k in BB){
	freq <- append(freq,k) 
}
}
S3 <- data.frame(site,aa,freq)
#S3$site <- gsub("S_", "", S3$site, fixed = T )#

S4 <- S3[!S3$aa == "x",]
S4$aa <- toupper(S4$aa)
S5 <- data.frame(Position=paste0(S4[,1],S4[,2]),S4)
S5 <- S5[order(S5$site,S5$freq),]
S5 <- data.frame(S5,lin=rep(linaje,dim(S5)[1]))

library(ggplot2)
ff <- ggplot(S5, aes(x=reorder(S5$Position,as.numeric(S5$site)), y=S5$freq, fill=site))
ff2 <- ff + geom_bar(stat="identity", color="black",size = 0.8) + theme_minimal() + scale_x_discrete(guide = guide_axis(angle = 45)) +
xlab("Position") + ylab("Frequency") + labs(title=paste(label,"mutations","run:",run,linaje)) + 
theme(legend.position = "none") + theme(axis.text=element_text(size=10)) + 
geom_text(aes(label=freq), vjust=-0.5, color="black", position = position_dodge(1), size=3.5)
ff2
print(S3)
print(S5)
print(ff2)
print(mutates_sites) 
x <- structure(3,class=c(S3,S5,ff2,mutates_sites))
return(x)
}

#########
### 4 ###
#########

#funcion para obtener plots de sitios variables por linaje#
#dichos sitios estan definidos si la variacion interna esta#
#presente en un numero minimo de genomas#

######## the function #######
#data = tab#
#linaje = "P.1.12"#
#genomas = 1#
#run = 31#
barras2 <- function(data,linaje,genomas,run,label,inic){
data = data #archivo fasta leido con read.fasta previamente#
linaje = linaje #linaje sobre el que se quiere realizar la estimacion#
run = run 
inic = inic
genomas = genomas
label = label
names = unique(data$lineage)

#modifica los siguiente (retira los # la proxima vez)#

#ch <- data[data$lineage == linaje,]#
ch <- data
ch2 <- ch[,inic:(ncol(data)-1)]

#### WE WILL ESTIMATE THE VARIABLE SITES ####

ngenomes <- genomas
#prepare the loop#
freq <- 0
r <- matrix(nrow=dim(ch2)[2],ncol=1)
table <- matrix(nrow=dim(ch2)[2],ncol=1)
for(i in 1:dim(ch2)[2]) {
freq[i] <- as.data.frame(summary(as.factor(ch2[,i])))
table[i] <- data.frame(n=rep(i,length(unlist(freq[[i]]))),changes=unlist(freq[[i]]))
r[i] <- as.data.frame(names(summary(as.factor(ch2[,i]))))
print(freq[[i]])
print(r[[i]])
}
freq
table
r

##### WE WILL MERGE THE PRODUCED DATA #####
s <- data.frame(position=unlist(table),aa=unlist(r),freq=unlist(freq))
s
names(s)
s2 <- s[!s$aa=="x" & s$freq >= ngenomes,]
s2[order(-s2$freq),]
names(s2)
#write.csv(s, file="frequencies_S.csv",row.names=FALSE)#

##### WE WILL ESTIMATE THE FREQUENCIES FOR VARIABLE SITES WITH MORE THAN 2 STATES #####
s3 <- data.frame(s2, n=rep(1,dim(s2)[1]))
s3[1:25,]
freq <- aggregate(s3$n, by=list(s3$position), FUN=sum)
names(freq) <- c("position","freq")
freq2 <- freq[order(-freq$freq),]
names(freq2)
lists <- freq[freq$freq >= 2,]
lists[,1]
dim(lists)
length(lists[,1])
#write.csv(lists, file=paste0("sites_in",genomas,"genomes",".csv"),row.names=FALSE)#

S <- ch[order(ch$taxon),c(1:(inic-1),lists[,1]+(inic-1))]
names(S)
dim(S)
mutates_sites <- names(S)
write.csv(S, file=paste0(label,"_selected_positions.csv"), row.names=FALSE)

len <- length(lists[,1])

S2 <- data.frame(S[,c(inic:dim(S)[2])],n=rep(1,dim(S)[1]))

names(S2)
print(S2)
v <- names(S[inic])
nam <- function(len){
if(len == 1) return(c(v,"n"))
else if(len >= 2) return(names(S2))
}

b <- nam(len)
names(S2) <- b
S2

site <- c()
aa <- c()
freq <- c()

for (i in 1:(ncol(S2)-1)){
	AA <- aggregate(S2$n, by=list(S2[,i]), FUN=sum)[,1]
	BB <- aggregate(S2$n, by=list(S2[,i]), FUN=sum)[,2]
	site <- append(site,rep(colnames(S2)[i],length(AA)))
for (j in AA){
	aa <- append(aa, j)
}
for(k in BB){
	freq <- append(freq,k) 
}
}
S3 <- data.frame(site,aa,freq)
#S3$site <- gsub("S_", "", S3$site, fixed = T )#

S4 <- S3[!S3$aa == "x",]
S4$aa <- toupper(S4$aa)
S5 <- data.frame(Position=paste0(S4[,1],S4[,2]),S4)
S5 <- S5[order(S5$site,S5$freq),]
S5 <- data.frame(S5,lin=rep(linaje,dim(S5)[1]))
write.csv(S5, file=paste0(label,"_positions.csv"), row.names=FALSE)

library(ggplot2)
ff <- ggplot(S5, aes(x=reorder(S5$Position,as.numeric(S5$site)), y=S5$freq, fill=site))
ff2 <- ff + geom_bar(stat="identity", color="black",size = 0.8) + theme_minimal() + scale_x_discrete(guide = guide_axis(angle = 45)) +
xlab("Position") + ylab("Frequency") + labs(title=paste(label,"mutations","run:",run,linaje)) + 
theme(legend.position = "none") + theme(axis.text=element_text(size=10)) + 
geom_text(aes(label=freq), vjust=-0.5, color="black", position = position_dodge(1), size=3.5)
ff2
print(S3)
print(S5)
print(ff2)
print(mutates_sites) 
x <- structure(3,class=c(S3,S5,ff2,mutates_sites))
return(x)
write.csv(S5, file=paste0(label,"positions.csv"), row.names=FALSE)
}
