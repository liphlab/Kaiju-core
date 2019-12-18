#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
#°°°°°°°°°°°°°°°°°°°°°°°°°° Mock 1 °°°°°°°°°°°°°°°°°°°°°°°°°#
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#

Kaiju<-read.delim("Mock1_ALLgenus.txt",na.strings="-",header=FALSE)
KaijuReads<-as.numeric(as.character(Kaiju[3:(nrow(Kaiju)-5),2]))
KaijuRank<-as.character(Kaiju[3:(nrow(Kaiju)-5),3])
KaijuTax<-matrix(0,ncol=6,nrow=length(KaijuReads))
for (i in 1:length(KaijuReads))
{KaijuTax[i,2:6]<-unlist(strsplit(KaijuRank[i],split=";"))}
KaijuTax[,1]<-rep("Bacteria",length(KaijuReads))
KaijuMat<-cbind(KaijuReads,KaijuTax)

Na.contained<-unique(sort(which(KaijuMat[,]=="NA", arr.ind=TRUE)[,1]))
KaijuMat<-KaijuMat[-Na.contained,]

rm(Kaiju,KaijuRank,KaijuReads,KaijuTax,i,Na.contained)


Kaiju<-read.delim("Mock1_PFAMgenus.txt",na.strings="-",header=FALSE)
KaijuReads<-as.numeric(as.character(Kaiju[3:(nrow(Kaiju)-5),2]))
KaijuRank<-as.character(Kaiju[3:(nrow(Kaiju)-5),3])
KaijuTax<-matrix(0,ncol=6,nrow=length(KaijuReads))
for (i in 1:length(KaijuReads))
{KaijuTax[i,2:6]<-unlist(strsplit(KaijuRank[i],split=";"))}
KaijuTax[,1]<-rep("Bacteria",length(KaijuReads))
KaijuMatPFAM<-cbind(KaijuReads,KaijuTax)

Na.contained<-unique(sort(which(KaijuMatPFAM[,]=="NA", arr.ind=TRUE)[,1]))
KaijuMatPFAM<-KaijuMatPFAM[-Na.contained,]

rm(Kaiju,KaijuRank,KaijuReads,KaijuTax,i,Na.contained)


GenusListTrue<-c("BifidobacteriaceaeBifidobacterium",
"BacteroidaceaeBacteroides","ClostridiaceaeClostridium",
"LactobacillaceaeLactobacillus","EnterobacteriaceaeEscherichia",
"StreptococcaceaeStreptococcus","EnterobacteriaceaeSalmonella")
TrueReads<-c(15.88,18.36,9.53,24.75,6.42,8.32,16.74)

GenusListKaiju<-rep("0",nrow(KaijuMat))
for (i in 1:nrow(KaijuMat))
GenusListKaiju[i]<-paste(KaijuMat[i,c(6,7)],collapse="")

GenusListPFAM<-rep("0",nrow(KaijuMatPFAM))
for (i in 1:nrow(KaijuMatPFAM))
GenusListPFAM[i]<-paste(KaijuMatPFAM[i,c(6,7)],collapse="")

GenusList<-unique(sort(c(GenusListKaiju,GenusListPFAM)))
GenusMat<-matrix(0,nrow=length(GenusList),ncol=4)
colnames(GenusMat)<-c("FamilyGenus","True","Kaiju","PFAM-Kaiju")

GenusMat[,1]<-GenusList
for (i in 1:nrow(GenusMat))
{if(length(which(GenusListTrue==GenusList[i]))>0)
    GenusMat[i,2]<-TrueReads[which(GenusListTrue==GenusList[i])]
    if(length(which(GenusListKaiju==GenusList[i]))>0)
    GenusMat[i,3]<-sum(as.numeric(KaijuMat[which(GenusListKaiju==GenusList[i]),1]))
    if(length(which(GenusListPFAM==GenusList[i]))>0)
    GenusMat[i,4]<-sum(as.numeric(KaijuMatPFAM[which(GenusListPFAM==GenusList[i]),1]))
}

GenusMat<-GenusMat[order(GenusMat[,1]),]

rm(GenusListKaiju,GenusListPFAM,GenusListTrue,GenusList,i,TrueReads,KaijuMat,KaijuMatPFAM)

AbundanceMat<-GenusMat[,2:4]
class(AbundanceMat)<-"numeric"
row.names(AbundanceMat)<-GenusMat[,1]
rm(GenusMat)

# If you wish to compare more taxonomic methods, add columns to the previous
# matrix in a similar manner so to obtain a Gx(1+M) matrix where rows are genera
# and columns are methods. In our case M=9 (True, Kaiju, Core-Kaiju, MetaPhlAn2,
# Dada 2 + Silva, Dada 2 + GG, Dada 2 + RDP, Qiime2 + SILVA, Qiime 2 + GG). The
# additional first column is neede to label genera (rows' names).

# Check whether same genera are labelled with different names via different
#  methods, in which case you need to sum up the abundances
IndexPept<-which(row.names(AbundanceMat)=="PeptostreptococcaceaePaeniclostridium")
IndexClo<-which(row.names(AbundanceMat)=="ClostridiaceaeClostridium")
AbundanceMat[IndexClo,]<-AbundanceMat[IndexClo,]+AbundanceMat[IndexPept,]
AbundanceMat<-AbundanceMat[-IndexPept,]

IndexPept<-which(row.names(AbundanceMat)=="PeptostreptococcaceaeParaclostridium")
IndexClo<-which(row.names(AbundanceMat)=="ClostridiaceaeClostridium")
AbundanceMat[IndexClo,]<-AbundanceMat[IndexClo,]+AbundanceMat[IndexPept,]
AbundanceMat<-AbundanceMat[-IndexPept,]

IndexPept<-which(row.names(AbundanceMat)=="Peptostreptococcaceae[Clostridium]")
IndexClo<-which(row.names(AbundanceMat)=="ClostridiaceaeClostridium")
AbundanceMat[IndexClo,2:ncol(AbundanceMat)]<-AbundanceMat[IndexClo,2:ncol(AbundanceMat)]+AbundanceMat[IndexPept,2:ncol(AbundanceMat)]
AbundanceMat<-AbundanceMat[-IndexPept,]

IndexPept<-which(row.names(AbundanceMat)=="PeptostreptococcaceaeClostridium_XI")
IndexClo<-which(row.names(AbundanceMat)=="ClostridiaceaeClostridium")
AbundanceMat[IndexClo,2:ncol(AbundanceMat)]<-AbundanceMat[IndexClo,2:ncol(AbundanceMat)]+AbundanceMat[IndexPept,2:ncol(AbundanceMat)]
AbundanceMat<-AbundanceMat[-IndexPept,]

IndexPept<-which(row.names(AbundanceMat)=="Clostridiaceae Clostridiumsensustricto")
IndexClo<-which(row.names(AbundanceMat)=="ClostridiaceaeClostridium")
AbundanceMat[IndexClo,2:ncol(AbundanceMat)]<-AbundanceMat[IndexClo,2:ncol(AbundanceMat)]+AbundanceMat[IndexPept,2:ncol(AbundanceMat)]
AbundanceMat<-AbundanceMat[-IndexPept,]

IndexPept<-which(row.names(AbundanceMat)=="PeptostreptococcaceaeClostridium")
IndexClo<-which(row.names(AbundanceMat)=="ClostridiaceaeClostridium")
AbundanceMat[IndexClo,2:ncol(AbundanceMat)]<-AbundanceMat[IndexClo,2:ncol(AbundanceMat)]+AbundanceMat[IndexPept,2:ncol(AbundanceMat)]
AbundanceMat<-AbundanceMat[-IndexPept,]


# Apply PFAM protocol: retain genera with abundance of at least 11 reads in PFAM-Kaiju column
EstVsTr<-AbundanceMat[,c(1,2)]
PFAMSelection<-which(AbundanceMat[,3]>10)
TrueSelection<-which(AbundanceMat[,1]>0)
UnionSelection<-unique(sort(c(PFAMSelection,TrueSelection)))
EstVsTr<-EstVsTr[UnionSelection,]

# Get relative abundance vector (RSA)
EstVsTr[,1]<-EstVsTr[,1]/sum(EstVsTr[,1])
EstVsTr[,2]<-EstVsTr[,2]/sum(EstVsTr[,2])

# Compare via a linear fit the two RSAs (true and predicted)
summary(lm(EstVsTr[,2]~EstVsTr[,1] ))

