#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
#°°°°°°°°°°°°°° CAMI high-complexity sample 1 °°°°°°°°°°°°°°#
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°#
# Load the true genus abundances provided at the CAMI site
gs<-read.table("goldstandard_high_1_genus.txt", sep="\t", header = T, stringsAsFactors = F)
gs<-as.numeric(gs[,5])

TrueMat<-matrix("0",nrow=length(gs), ncol=5)
for (i in 1:length(gs))
{TrueMat[i,]<-unlist(strsplit(TaxaLab[i], split='|', fixed=T))[2:6]}

ListTrueTaxa<-c()
for (i in 1:nrow(TrueMat))
ListTrueTaxa<-c(ListTrueTaxa,TrueMat[i,5])

# Load results for Kaiju 1.0
Kaiju<-read.delim("CH1.txt",na.strings="-",header=FALSE)
KaijuReads<-as.numeric(as.character(Kaiju[3:(nrow(Kaiju)-5),2]))
KaijuRank<-as.character(Kaiju[3:(nrow(Kaiju)-5),3])
KaijuTax<-matrix(0,ncol=6,nrow=length(KaijuReads))
for (i in 1:length(KaijuReads))
{KaijuTax[i,2:6]<-unlist(strsplit(KaijuRank[i],split=";"))}
KaijuTax[,1]<-rep("Bacteria",length(KaijuReads))
KaijuMat<-cbind(KaijuReads,KaijuTax)

Na.contained<-unique(sort(which(KaijuMat[,]=="NA", arr.ind=TRUE)[,1]))
if (length(Na.containes)>0) KaijuMat<-KaijuMat[-Na.contained,]

rm(Kaiju,KaijuRank,KaijuReads,KaijuTax,i,Na.contained)

ListTaxaKaiju<-c()
for (i in 1:nrow(KaijuMat))
ListTaxaKaiju<-c(ListTaxaKaiju,KaijuMat[i,7])

# Load results for Kaiju against PFAM database
Kaiju<-read.delim("CH1_PFAM.txt",na.strings="-",header=FALSE)
KaijuReads<-as.numeric(as.character(Kaiju[3:(nrow(Kaiju)-5),2]))
KaijuRank<-as.character(Kaiju[3:(nrow(Kaiju)-5),3])
KaijuTax<-matrix(0,ncol=6,nrow=length(KaijuReads))
for (i in 1:length(KaijuReads))
{KaijuTax[i,2:6]<-unlist(strsplit(KaijuRank[i],split=";"))}
KaijuTax[,1]<-rep("Bacteria",length(KaijuReads))
KaijuMatPFAM<-cbind(KaijuReads,KaijuTax)

Na.contained<-unique(sort(which(KaijuMatPFAM[,]=="NA", arr.ind=TRUE)[,1]))
if (length(Na.containes)>0) KaijuMatPFAM<-KaijuMatPFAM[-Na.contained,]

rm(Kaiju,KaijuRank,KaijuReads,KaijuTax,i,Na.contained)

# Obtain Core-Kaiju results by applying the absolute threshold
KaijuPFAMAbb<-as.numeric(KaijuMatPFAM[,1])
KaijuMatPFAM<-KaijuMatPFAM[which(KaijuPFAMAbb>20),]

ListTaxaKaijuPFAM<-c()
for (i in 1:nrow(KaijuMatPFAM))
ListTaxaKaijuPFAM<-c(ListTaxaKaijuPFAM,KaijuMatPFAM[i,7])

PFAMprot<-which(ListTaxaKaiju%in%ListTaxaKaijuPFAM)

CoreKaijuMat<-KaijuMat[PFAMprot,]
rm(KaijuPFAMAbb, KaijuMatPFAM, ListTaxaKaijuPFAM, PFAMprot)


# Put all found genera together
ListTaxa<-unique(sort(c(ListTrueTaxa, ListTaxaKaiju)))
rm(ListTaxaKaiju, ListTrueTaxa)

TaxaMat<-matrix(0,ncol=3, nrow=length(ListTaxa))
colnames(TaxaMat)<-c("True","Core-Kaiju","Kaiju")
rownames(TaxaMat)<-ListTaxa

for (i in 1:nrow(TrueMat))
{taxa_i<-TrueMat[i,5]
    index_taxa<-which(ListTaxa==taxa_i)
    TaxaMat[index_taxa,1]<-TaxaMat[index_taxa,1]+gs[i]}

for (i in 1:nrow(CoreKaijuMat))
{taxa_i<-CoreKaijuMat[i,7]
    index_taxa<-which(ListTaxa==taxa_i)
    TaxaMat[index_taxa,2]<-TaxaMat[index_taxa,2]+as.numeric(CoreKaijuMat[i,1])}

for (i in 1:nrow(KaijuMat))
{taxa_i<-KaijuMat[i,7]
    index_taxa<-which(ListTaxa==taxa_i)
    TaxaMat[index_taxa,3]<-TaxaMat[index_taxa,3]+as.numeric(KaijuMat[i,1])}

TaxaMat[,1]<-TaxaMat[,1]/sum(TaxaMat[,1])
TaxaMat[,2]<-TaxaMat[,2]/sum(TaxaMat[,2])
TaxaMat[,3]<-TaxaMat[,3]/sum(TaxaMat[,3])

rm(CoreKaijuMat,gs,KaijuMat,TrueMat,i, index_taxa, ListTaxa, taxa_i)


###################### Core-Kaiju results #####################
TruePos<-length(which(TaxaMat[,1]!=0 & TaxaMat[,2]!=0))
FalsePos<-length(which(TaxaMat[,1]==0 & TaxaMat[,2]!=0))
FalseNeg<-length(which(TaxaMat[,1]!=0 & TaxaMat[,2]==0))

length(which(TaxaMat[,2]!=0)) # prediction for G

# Precision (true pos/ true pos+ false pos)
Prec<-TruePos/(TruePos+FalsePos)

# Recall (true pos/true pos+false neg)
Rec<-TruePos/(TruePos+FalseNeg)

# F1 Score 2*prec*recall/(prec+recall)
F1score<-2*Prec*Rec/(Prec+Rec)

# Linear fit between true and predicted abundances
EstVsTr<-TaxaMat[,c(1,2)]
EstVsTr<-EstVsTr[-which(EstVsTr[,1]==0& EstVsTr[,2]==0),]
summary(lm(EstVsTr[,2]~EstVsTr[,1] ))

plot(EstVsTr, pch=19, xlab="True", col="red",ylab="Core-Kaiju")
abline(lm(EstVsTr[,2]~EstVsTr[,1] ), col="green", lwd=2)
abline(a=0,b=1, col="black", lwd=2,lty=3)
