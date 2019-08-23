library(dplyr)
library(ggplot2)
library("beeswarm")  
library("rms")  # ORM ######
###### data1: APOBEC deletion data from APOBEC NG paper  ##########

ng.del<-read.table("ng.TCGA.del")
ng.del$V2<-substr(ng.del$V2,1,12)
ng.del$V2<-gsub("-",".",ng.del$V2,perl=TRUE)
ng.del<-ng.del[order(ng.del$V2),]
ng.del<-ng.del[!duplicated(ng.del$V2),]
row.names(ng.del)<-ng.del$V2

x.t1<-ng.del[ng.del$V8 <2,]

x.t1$V8[x.t1$V8 ==1] <- -1
x.t1$V8[x.t1$V8 ==0] <- -2
x.t1$V8[x.t1$V8 ==2] <- 0

##### APOBEC signature data ######

dat<-read.table("pancancer.maf_sorted_sum_all_fisher_Pcorr.txt",head=TRUE,sep="\t")
dat.sub<-dat[,c(1,6,53)] 
dat.sub$Sample<-substr(dat.sub$Sample,1,12)
dat.sub$Sample<-gsub("-",".",dat.sub$Sample,perl=TRUE)
dat.sub<-dat.sub[order(dat.sub$Sample),]
dat.sub<-dat.sub[!duplicated(dat.sub$Sample),]
######## end 1 #####################

#names(clin)<-toupper(names(clin))
#rownames(clin)<-clin$BCR_PATIENT_BARCODE
#clin<-clin[,-1]
#clin<-t(clin)

#clin<-as.data.frame(clin)

######## data2: clin -> race  #################
#clin<-read.table("race.list",sep="\t",stringsAsFactors=FALSE)
#clin<-read.table("race.list",sep="\t")
#names(clin)<-c("sample","race")
#clin<-clin[order(clin$sample),]
#clin$sample<-gsub("-",".",clin$sample,perl=TRUE)
################## end 2 #######



###### data3:  APOBEC3 CNA from TCGA ######
cna<-read.table("APOBEC.data_CNA.txt",head=TRUE)

#cna.apo<-cna[,c(2,23517,23518)]

#cna.apo.rm<-cna.apo[apply(cna.apo[,-1],1,sum) == 0,]

cna$Patient<-substr(cna$Patient,1,12) 
cna$Patient<-gsub("-",".",cna$Patient,perl=TRUE) 
cna<-cna[order(cna$Patient),]
names(cna)<-c("Patient","sAPOBEC3A","sAPOBEC3B")
###### end 1  ##########

####### data4: APOBEC3 germline CN from TCGA #######3

gcna1<-read.table("pancancer.seg.input1",head=FALSE)
names(gcna1)<-c("sample","chr","s1","s2","probes","seg")
#cna.apo<-cna[,c(2,23517,23518)]

gcna1$sample<-substr(gcna1$sample,1,12) 
gcna1$sample<-gsub("-",".",gcna1$sample,perl=TRUE) 

gcna1<- gcna1[gcna1$probes >50,]
gcna1<-aggregate(seg~sample,data=gcna1,median) 
gcna1<-gcna1[gcna1$seg > -0.2,]
#cna.apo.rm<-cna.apo[apply(cna.apo[,-1],1,sum) == 0,]

gcna<-read.table("pancancer.seg.input",head=FALSE)
names(gcna)<-c("sample","chr","s1","s2","probes","seg")

#cna.apo<-cna[,c(2,23517,23518)]

gcna$sample<-substr(gcna$sample,1,12) 
gcna$sample<-gsub("-",".",gcna$sample,perl=TRUE) 
#gcna<- gcna[gcna$seg < -0.5  & gcna$probes < 50,]
#gcna1<- gcna[gcna$seg >  -0.1,]
gcna2<-aggregate(seg~sample,data=gcna,median)

#gcna2<- rbind(gcna1[gcna1$probes > 50,],gcna)

#cna.apo.rm<-cna.apo[apply(cna.apo[,-1],1,sum) == 0,]


####### end 4 ###################

####### data5: APOBEC expression data from TCGA #######
#expr<-read.table("APOBEC.iso.expr",head=TRUE,row.names=1,stringsAsFactors=FALSE)
#expr<-t(expr)
#write.table(expr,"APOBEC.iso.expr",quote=F,sep="\t")
Ctype<-read.table("Pancancer_expr_NM_APOBEC.txt",head=TRUE,stringsAsFactors=FALSE)

Ctype$Samples<-substr(Ctype$Samples,1,12) 
Ctype<-Ctype[!duplicated(Ctype$Samples),]

expr<-read.table("APOBEC.iso.expr",head=TRUE,row.names=1,stringsAsFactors=FALSE)
x1<-rownames(expr)
expr1<-cbind(x1,expr)
expr1$x1<-substr(expr1$x1,1,12) 
expr1<-expr1[!duplicated(expr1$x1),]
rownames(expr1)<-expr1$x1
expr<-expr1[,-1]
#expr<-log2(expr1[,-1]+0.0001)
###### end data5 ############

#### data6: nenotigene from TCIA  ###########
neo<-read.table("TCIA/neoantigensAll.tsv",head=TRUE,sep="\t")

neo.acc<-data.frame(table(neo$patientBarcode))
neo.acc$Var1<-gsub("-",".",neo.acc$Var1,perl=TRUE) 
neo.acc<-neo.acc[order(neo.acc$Var1),]

##### main analysis #########
##### anaysis 1: construct matrix ######

fx<- 1
Ctypes<-c("BLCA", "BRCA", "CESC", "LUAD", "LUSC", "HNSC", "STAD","PAAD","THCA", "KIRP")
#Ctypes<-c("BLCA", "BRCA", "CESC", "LUAD", "LUSC", "HNSC", "THCA", "KIRP")
######## section I: Investigate expression level with signature ########
com<-intersect(dat.sub$Sample,intersect(Ctype$Samples,rownames(expr)))
com<-intersect(dat.sub$Sample,intersect(Ctype$Samples,rownames(expr)))

Ctype<-Ctype[Ctype$Samples %in% com,]
Ctype<-Ctype[order(Ctype$Samples),]
Ctype<-Ctype[,-1]
dat.sub<-dat.sub[dat.sub$Sample %in% com,]
dat.sub<-dat.sub[!duplicated(dat.sub$Sample),]
expr<-expr[rownames(expr) %in% com,]
expr<-expr[order(rownames(expr)),]


mat<-cbind(dat.sub,Ctype,expr)
mat<-mat[mat$Cancer %in% Ctypes,]

Sig<-log2(1+mat$mutations * mat$X.tCw_to_G.tCw_to_T._per_mut)
mat<-cbind(mat,Sig)
rownames(mat)<-mat$Sample
### all and subtype analysis 1: ading pam50 ###
pam50<-read.table("BRCA_pam50.txt")
pam50$V1 <- substr(pam50$V1,1,12)
pam50$V1 <- gsub("-",".",pam50$V1,perl =TRUE)
# statistical analysis #######
for( type in Ctypes[1:10])
{
print(type) 
df = subset(mat[mat$Cancer %in% paste(type),],sel=c(12:17,27)); colnames(df) = c("x1","x2","x3","x4","x5","x6","y")
df = subset(mat[mat$Cancer %in% paste(type),],sel=c(12:17,3)); colnames(df) = c("x1","x2","x3","x4","x5","x6","y")
df = subset(mat[mat$Cancer %in% paste(type),],sel=c(24:27)); colnames(df) = c("x1","x2","x3","y")
df = subset(mat[mat$Cancer %in% paste(type),],sel=c(5,6,27)); colnames(df) = c("x1","x2","y")
df$x1 <- log2(df$x1+fx)
df$x2 <- log2(df$x2+fx)
df$x3 <- log2(df$x3+fx)
df$x4 <- log2(df$x4+fx)
df$x5 <- log2(df$x5+fx)
df$x6 <- log2(df$x6+fx)

#### subtype #####
comSam<- intersect(rownames(df),pam50$V1)
pam50<- pam50[!duplicated(pam50$V1),]
rownames(pam50) <- pam50$V1
pam50<- pam50[pam50$V1 %in% comSam,]
pam50run <- pam50[rownames(df),]$V2
df1<-cbind(df,pam50run)
#df$y <- log2(df$y)
pam<-c("Basal","Her2","LumA","LumB") 
x<-orm( y ~ x1+x2+x3 + x4 + x5 + x6,data=df1[df1$pam50run%in% "Her2",],family=probit)
x<-orm( y ~ x1,data=df1[df1$pam50run%in% "Her2",],family=probit)
x<-orm( y ~ x1+x2+x3 + x4 + x5 + x6,data=df1[df1$pam50run%in% paste(pam[2]),],family=probit)
x<-orm( y ~ x1+x2+x3 + x4 + x5 + x6,data=df1[df1$pam50run%in% paste(pam[3]),],family=probit)
x<-orm( y ~ x1+x2+x3 + x4 + x5 + x6,data=df1[df1$pam50run%in% paste(pam[4]),],family=probit)
x<-orm( y ~ x1+x2+x3 + x4 + x5 + x6,data=df,family=probit)
x<-orm( y ~ x1+x2+x3 + x4 + x5 + x6,data=df,family=probit)

x<-orm( y ~ x4,data=df1[df1$pam50run%in% paste(pam[4]),],family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]   
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2]))) 

pdf(paste(type,"isoform.expr.pdf",sep="."))
boxplot(df[1:6])
dev.off()

#print(summary(lm( y ~ x1 , data=df))$coefficients[2,c(1,4)])
x<-orm( y ~ x1,data=df,family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]   
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2]))) 
x<-orm( y ~ x2,data=df,family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]   
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2]))) 
x<-orm( y ~ x3,data=df,family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]   
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2]))) 
x<-orm( y ~ x4,data=df,family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]   
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2]))) 
x<-orm( y ~ x5,data=df,family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]   
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2]))) 
x<-orm( y ~ x6,data=df,family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]   
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2]))) 
}

######### end section I ######
######## section II: Investigate neoantigen with signature ########
#com<-intersect(dat.sub$Sample,intersect(Ctype$Samples,neo.acc$Var1))
com<-intersect(mat$Sample,neo.acc$Var1)
mat.1<-cbind(mat[mat$Sample %in% com,],neo.acc[neo.acc$Var1 %in% com,])


for( type in Ctypes)
{
print(type) 

df = subset(mat.1[mat.1$Cancer %in% paste(type),],sel=c(27,29)); colnames(df) = c("x","y")

print(summary(lm( y ~ x , data=df))$coefficients[2,c(1,4)])
}
########## end II #############




gcna2<-gcna2[!duplicated(gcna2$sample),]
gcna2<-gcna2[order(gcna2$sample),]
#





##### directly start from here ########3
rownames(mat.1)<-mat.1$Sample
mat<-mat.1
x<-as.data.frame(cbind(rownames(mat),5))
names(x)<-c("sample","value")
x$sample<-as.character(x$sample)

gcna2<-gcna2[gcna2$sample %in% mat$Sample,]
x$value<-as.numeric(as.character(x$value))   
x$value[x$sample %in% gcna2$sample] <- gcna2$seg

#x.t<-x.t[order(x.t$V2),]
#x.t<-x.t[x.t$V2%in% x$sample,]


x.t<-x[,-1]
mat.n<-cbind(mat,x.t)
#mat.n$x.t<-as.numeric(as.character(mat.n$x.t))
mat.n$x.t[mat.n$x.t < -1] <- -2
mat.n$x.t[mat.n$x.t < -0.2 & mat.n$x.t >= -1 ] <- -1
mat.n$x.t[mat.n$x.t > -0.2 & mat.n$x.t <= 5 ] <- 0
#mat.n<-mat.n[rownames(mat.n) %in% gcna3$sample,] 
mat.n<-mat.n[mat.n$x.t <= 0,]

x.t1<-x.t1[x.t1$V2 %in% mat.n$Sample,]
mat.n$x.t[mat.n$Sample %in% x.t1$V2] <- as.numeric(x.t1$V8)
mat.n<-mat.n[rownames(mat.n) %in% c(gcna1$sample,gcna2$sample),] 


com<-intersect(cna$Patient,mat.n$Sample)


cna<-cna[cna$Patient %in% com,]
cna<-cna[!duplicated(cna$Patient),]
cna<-cna[,-1]

mat.n<-mat.n[mat.n$Sample%in% com,]

mat.n2<-cbind(mat.n,cna)


# statistical analysis #######3

for( type in Ctypes )
{
print(type) 
df = subset(mat.n[mat.n$Cancer %in% paste(type),],sel=c(18,22)); colnames(df) = c("x","y")
print(summary(lm( y ~ x, data=df))$coefficients[2,c(1,4)])

}
###### total mutations vs neoantigen ######
pdf("Neoantigen_signature.pdf")
#X[X$Freq > 100,]$Var
for( type in Ctypes )
{
print(type)
df = subset(mat[mat$Cancer %in% paste(type),],sel=c(3,13)); colnames(df) = c("x","y")

#df$y<-log2df$y)
df$y<-log2(df$y)
print(ggplot(data=df,aes(x,y)) + 
  stat_density2d(aes(fill= ..level..),geom='polygon', contour = TRUE) + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=glm,linetype=2,colour="gray",se=F) + 
  guides(alpha="none") +
  geom_point(alpha = 1/5))
} 
dev.off()


###### total mutations vs neoantigen | consider germline/somatic ######
for(i in c(12,14,15))
{
pdf(paste(paste("signature_ISO",i,sep=""),".pdf",sep=""))
#X<-as.data.frame(table(mat.n$Cancer))
#X[X$Freq > 100,]$Var
for( type in  Ctypes)
{
print(type)
df = subset(mat[mat$Cancer %in% paste(type),],sel=c(i,27))
colnames(df) = c("x","y")
#df<-df[df$x >0,]
#X2<-mat.n[mat.n$Cancer %in% paste(type) & mat.n$x.t == 0,]
#xa<-apply(X1[,1:3],1,sum)
#df<-cbind(xa,X

colnames(df) = c("x","y")
df<-as.data.frame(df)
df$x<-log2(df$x+fx)
#df$y<-log2(df$y)
#df$x<-log2(df$x)
print(ggplot(data=df,aes(x,y)) + 
  stat_density2d(aes(fill= ..level..),geom='polygon', contour = TRUE) + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=glm,linetype=2,colour="gray",se=F) + 
  guides(alpha="none") +
  geom_point(alpha = 1/5))
} 
dev.off()
}



#### isoforms vs deletion ######

hap<-read.table("APOBECH.hap1", head=TRUE)
hap$Sample<- gsub("-",".",hap$Sample,perl = TRUE)
pdf("Updated1.load_deletion.BC.pdf")
for( type in Ctypes)
{
print(type)

Xt<-mat.n2[mat.n2$Cancer %in% paste(type) & mat.n2$sAPOBEC3A <=0,]
mA3B<-median(Xt$APOBEC3B)

Xt<-Xt[(Xt$x.t < 0 & Xt$APOBEC3B < mA3B) | Xt$x.t ==0,]
Xtdel <- Xt[Xt$x.t < 0,]

write.table(Xtdel[,c(1,14)],paste(type,"del",sep="."),quote=F,sep="\t")
}


print(summary(lm(log2(Freq) ~ x.t,data=Xt))$coeff[2,c(1,4)])
print(summary(lm(log2(fx+uc003awn) ~ x.t,data=Xt))$coeff[2,c(1,4)])
print(summary(lm(log2(fx+uc011aoc) ~ x.t,data=Xt))$coeff[2,c(1,4)])
print(summary(lm(log2(fx+uc003awo) ~ x.t,data=Xt))$coeff[2,c(1,4)])


###### #####

x<-orm( X.tCw_to_G.tCw_to_T._per_mut  ~ x.t + tmp,data=Xt,family=probit)
x<-orm( X.tCw_to_G.tCw_to_T._per_mut  ~ tmp,data=Xt,family=probit)
x<-orm( Sig  ~ x.t+tmp ,data=Xt,family=probit)
x<-orm( Sig  ~ tmp,data=Xt,family=probit)
x<-orm( X.tCw_to_G.tCw_to_T._per_mut  ~ x.t,data = Xt,family=probit)
#x<-orm(X.tCw_to_G.tCw_to_T._per_mut  ~ x.t+log2(fx+uc011aoc)+x.t:log2(fx+uc011aoc),data=Xt,family=probit)
i1<-length(x$coeff); i2 <- i1-1
a<- x$coeff[i2:i1]
b<- sqrt(diag(x$var))
print(a[2])
print(2-2*pnorm(abs(a[2]/b[2])))

#boxplot( log2(Freq)  ~ x.t, data = Xt , add = T)
}
dev.off()


comI.name<-names(com1)
for(i in 34:55)
{
pdf(paste(paste("Immune_del",i,sep=""),".pdf",sep=""))
#X<-as.data.frame(table(mat.n$Cancer))
#X[X$Freq > 100,]$Var
##### end isoforms ######
for( type in as.character(X$Var ))
{
print(type)
#beeswarm( X.tCw_to_G.tCw_to_T._per_mut   ~ x.t, data = mat.n[mat.n$Cancer %in% paste(type)& mat.n$sAPOBEC3A <=0 & ((mat.n$uc011aoc == 0 & mat.n$x.t ==0) | mat.n$x.t < 0) ,],
beeswarm( Sig  ~ x.t, data = mat.n[mat.n$Cancer %in% paste(type)& mat.n$sAPOBEC3A <=0,],
method = 'swarm',
#pch = 16, cex= 0.5, col=c("#E69F00", "#56B4E9", "#009E73"),
pch = 16, cex= 0.5, col=c("#CC6666", "#9999CC", "#66CC99"),
xlab = '', ylab = '')
 
boxplot( log2(Freq)  ~ x.t, data = mat.n[mat.n$Cancer %in% paste(type) & mat.n$sAPOBEC3A <=0 ,], add = T)
}
dev.off()
}



########### Immune analysis after rm somatic CN APOBEC3  ######

ImmuneP<-read.table("Pancancer.Immune.cell.txt",head=TRUE,,stringsAsFactors=FALSE,sep="\t")

ImmuneP$Input.Sample<-substr(ImmuneP$Input.Sample,1,12)   
#ImmuneP<-ImmuneP[ImmuneP$Input.Sample %in% com,]
ImmuneP<-ImmuneP[!duplicated(ImmuneP$Input.Sample),]
ImmuneP<-ImmuneP[order(ImmuneP$Input.Sample),]
ImmuneP<-ImmuneP[ImmuneP$P.value < 0.2,]

comI<-intersect(mat.n2$Sample,ImmuneP$Input.Sample)

mat.In<-cbind(mat.n2[mat.n2$Sample %in% comI,],ImmuneP[ImmuneP$Input.Sample %in% comI,])



##### end isoforms ######
for(i in c(35:36,38,41,43))
{
pdf(paste(paste("Updated.Immune_del",i,sep=""),".pdf",sep=""))
for( type in Ctypes[6:10])
{
print(type)
type<-"BRCA"
Xt<-mat.In[mat.In$Cancer %in% paste(type),] #### for 2210 no deletion data
Xt<-mat.In[mat.In$Cancer %in% paste(type) & mat.In$sAPOBEC3A <=0,]
mA3B<-median(Xt$APOBEC3B)
Xt<-Xt[(Xt$x.t < 0 & Xt$APOBEC3B < mA3B) | Xt$x.t ==0,]
#beeswarm( X.tCw_to_G.tCw_to_T._per_mut   ~ x.t, data = mat.n[mat.n$Cancer %in% paste(type)& mat.n$sAPOBEC3A <=0 & ((mat.n$uc011aoc == 0 & mat.n$x.t ==0) | mat.n$x.t < 0) ,],
#beeswarm( comI.name[i]   ~ x.t, data = [mat.n$Cancer %in% paste(type)& mat.n$sAPOBEC3A <=0,],
df = subset(Xt,sel=c(30,i)); colnames(df) = c("x","y")
beeswarm(  y  ~ x, data = df,
method = 'swarm',
#pch = 16, cex= 0.5, col=c("#E69F00", "#56B4E9", "#009E73"),
pch = 16, cex= 0.5, col=c("#CC6666", "#9999CC", "#66CC99"),
xlab = '', ylab = '')
i=40
print(summary(lm( B.cells.naive ~ X.tCw_to_G.tCw_to_T._per_mut, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( B.cells.memory ~ X.tCw_to_G.tCw_to_T._per_mut, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD8 ~ X.tCw_to_G.tCw_to_T._per_mut, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD4.naive ~ X.tCw_to_G.tCw_to_T._per_mut, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD4.memory.activated ~ X.tCw_to_G.tCw_to_T._per_mut, data=Xt))$coefficients[2,c(1,4)])

print(summary(lm( B.cells.naive ~ uc011aoc, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( B.cells.memory ~ uc011aoc, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD8 ~ uc011aoc, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD4.naive ~ uc011aoc, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD4.memory.activated ~ uc011aoc, data=Xt))$coefficients[2,c(1,4)])

print(wilcox.test(Xt[Xt$x.t == 0, i],Xt[Xt$x.t < 0, i]))
print(summary(lm( B.cells.naive ~ x.t, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( B.cells.memory ~ x.t, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD8 ~ x.t, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD4.naive ~ x.t, data=Xt))$coefficients[2,c(1,4)])
print(summary(lm( T.cells.CD4.memory.activated ~ x.t, data=Xt))$coefficients[2,c(1,4)])
print(wilcox.test(Xt[Xt$x.t == 0, i],Xt[Xt$x.t < 0, i]))
boxplot( y  ~ x, data = df, add = T)
}
dev.off()
}







### for Immune vs signature, without considering deletion ###### ######
comI<-intersect(mat.1$Sample,ImmuneP$Input.Sample)

mat.In<-cbind(mat.1[mat.1$Sample %in% comI,],ImmuneP[ImmuneP$Input.Sample %in% comI,])
####### end #########
for(i in c(35:36,38,41,43))
{
i<-35
for( type in Ctypes)
{
print(type)
Xt<-mat.In[mat.In$Cancer %in% paste(type),]
xa<-log2(Xt[,29])
ya<-Xt[,i]
#print(orm( ya ~ xa, family = cloglog)$coefficients)
#print(orm( ya ~ xa, family = cloglog))
print(summary(lm( ya ~ xa))$coefficients[2,c(1,4)])  

}




###### total mutations vs neoantigen | consider germline/somatic ######
for(i in c(37,42))
{
pdf(paste(paste("Immune_Mutations.",i,sep=""),".pdf",sep=""))
pdf(paste(paste("Immune_signature.",i,sep=""),".pdf",sep=""))
#X<-as.data.frame(table(mat.n$Cancer))
#X[X$Freq > 100,]$Var
for( type in as.character(X$Var ))
{
print(type)
df = subset(mat.In[mat.In$Cancer %in% paste(type),],sel=c(22,i))
colnames(df) = c("x","y")
#df<-df[df$x >0,]
#X1<-mat.n[mat.n$Cancer %in% paste(type) & mat.n$x.t == 0,]
#xa<-apply(X1[,1:3],1,sum)
#df<-cbind(xa,X1[,18])
colnames(df) = c("x","y")
#df<-as.data.frame(df)
#df$y<-log10(df$y)
#df$x<-log10(df$x)
df$x<-log2(df$x)
print(ggplot(data=df,aes(x,y)) +
  stat_density2d(aes(fill= ..level..),geom='polygon', contour = TRUE) +
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=glm,linetype=2,colour="gray",se=F) +
  guides(alpha="none") +
  geom_point(alpha = 1/5))
}
dev.off()
}

