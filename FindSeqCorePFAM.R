# Upload data
list_proteomes<-read.table(file="list_proteomes_bacteria.txt",sep=" ",header=F,as.is=TRUE)
list_PFAM<-read.table(file="list_Pfam_bac.txt",sep=" ",header=F,as.is=TRUE)
load("matrix_CorePFAMvsPROTEOMES_bacteria")

# Core PFAM 
load(file="list_CorePFAM")
CorePFAM.names<-list_PFAM[CorePFAM,1]
n.CorePFAM<-length(CorePFAM)

# Create a list of n.CorePFAM elements, each consisting in a character
# matrix where we will store, for each proteomes containing it (rows), the names 
# of the proteins where the PFAM has been found (columns).
n.proteomes<-length(list_proteomes)

CorePROTEINS<-list()
for (i in 1:n.CorePFAM)
  CorePROTEINS[[i]]<-matrix(0,ncol=max(matrix_CorePFAM[i,]),nrow=n.proteomes)

CorePFAM.position<-list()
for (i in 1:n.CorePFAM)
  CorePFAM.position[[i]]<-matrix(0,ncol=max(matrix_CorePFAM[i,]),nrow=n.proteomes)

listPROTEINS<-c("0")

for (j in 1:n.proteomes)
  {prot.name<-list_proteomes[1,j]
   PROTfile<-read.table(file=paste("fold_bacteria/",prot.name,sep=""))

   for (i in 1:n.CorePFAM)
     if(matrix_CorePFAM[i,j]>0)
       {index_domain<-which(PROTfile[,6]==CorePFAM.names[i])
        CorePROTEINS[[i]][j,1:length(index_domain)]<-as.character(PROTfile[index_domain,1])
        CorePFAM.position[[i]][j,1:length(index_domain)]<-sapply(1:length(index_domain),function(x) paste(as.character(PROTfile[index_domain[x],c(2,3),]),collapse ="-"))
        listPROTEINS<-c(listPROTEINS,CorePROTEINS[[i]][j,])
        listPROTEINS<-unique(sort(listPROTEINS))
        listPROTEINS<-listPROTEINS[listPROTEINS!="0"]
     }
   
   if (j%%100==0) cat(j,", ",sep="")
   if (j%%1000==0) cat ("\n")
   
  }

# Check whether a Core PFAM can occurr more than once in the same protein
# and create a list of all proteines containing a Core PFAM
for (j in 1:n.proteomes)
  {for (i in 1:n.CorePFAM)
    if(matrix_CorePFAM[i,j]>1)
      {string<-CorePROTEINS[[i]][j,]
       string<-string[string!="0"]
       rep<-as.numeric(table(string))
       if(max(rep)>1) print(paste("Core PFAM",i, "occurrs",max(rep),"times in the same protein for proteome",j))}
  }

# Find nucleotidic sequences of PFAM using UniProt database
library(data.table)

# Create a list containing the nucleotidic sequence of each protein
# in listPROTEINS
n.proteins<-length(listPROTEINS)
listSEQ<-rep(0,n.proteins)

for (i in 1:n.proteins)
  {UniProtFile<-tryCatch(fread(paste("http://www.uniprot.org/uniprot/",listPROTEINS[i],".fasta",sep=""),sep="\t", showProgress = FALSE),error = function(e) {})
   if(length(UniProtFile)>0) 
     {UniProtSeq<-UniProtFile[1,]
   
   for (j in 2:nrow(UniProtFile))
     UniProtSeq<-paste(UniProtSeq,UniProtFile[j,],sep="")
   
   listSEQ[i]<-UniProtSeq}
  
   if (i%%100==0) cat(i,", ",sep="")
   if (i%%1000==0) cat ("\n")
}

# Check for missing proteins due to redundant proteome, thus eliminated
# from UniProt database. They are still available in UniParc. 
# (see https://www.uniprot.org/help/proteome_redundancy)
Error<-which(listSEQ=="0")

for (i in 1:length(Error))
{for (v in 1:100)
  {UniProtFile<-tryCatch(fread(
   paste("https://www.uniprot.org/uniprot/",listPROTEINS[Error[i]],".fasta?version=",v,sep=""),
   sep="\t", showProgress = FALSE),error = function(e) {})
   
   if(length(UniProtFile)==0) break 
  }
 
  if (v>1) 
    {UniProtFile<-tryCatch(fread(
     paste("https://www.uniprot.org/uniprot/",listPROTEINS[Error[i]],".fasta?version=",v-1,sep=""),
     sep="\t", showProgress = FALSE),error = function(e) {})}

  if(length(UniProtFile)>0) 
    {UniProtSeq<-UniProtFile[1,]

     for (j in 2:nrow(UniProtFile))
       UniProtSeq<-paste(UniProtSeq,UniProtFile[j,],collapse="")

    listSEQ[Error[i]]<-UniProtSeq}

   if (i%%10==0) cat(i,", ",sep="")
   if (i%%100==0) cat ("\n")
  }

# Create a list containing the nucleotidic sequence of each PFAM
# in each proteome
SeqPFAM<-list()
for (i in 1:n.CorePFAM)
  SeqPFAM[[i]]<-matrix(0,ncol=max(matrix_CorePFAM[i,]),nrow=n.proteomes)

for (j in 1:n.proteomes)
  {for (i in 1:n.CorePFAM)
    {Proteins<-unlist(lapply(CorePROTEINS[[i]][j,], function(x) {return(which(listPROTEINS==x))}))
     n.prot<-length(Proteins)
     
     if (n.prot>0)
      {Seq<-listSEQ[Proteins]
       Position<-as.numeric(unlist(strsplit(CorePFAM.position[[i]][j,][1:n.prot],split="-"))) 
       Beginning<-Position[seq(1,n.prot*2,2)]
       End<-Position[seq(2,n.prot*2,2)]
       SubSeq<-unlist(lapply(1:n.prot,function(x) {return(paste(unlist(strsplit(Seq[x],split=""))[Beginning[x]:End[x]],collapse=""))}))
       SeqPFAM[[i]][j,1:n.prot]<-SubSeq}
      }
    
    if (j%%100==0) cat(j,", ",sep="")
    if (j%%1000==0) cat ("\n")
  }

# Create a list containing the length of the nucleotidic sequence of
# each PFAM in each proteome
LengthSeqPFAM<-list()
for (i in 1:n.CorePFAM)
  LengthSeqPFAM[[i]]<-matrix(0,ncol=max(matrix_CorePFAM[i,]),nrow=n.proteomes)

for (j in 1:n.proteomes)
  {for (i in 1:n.CorePFAM)
    {String<-SeqPFAM[[i]][j,]
     LengthSeqPFAM[[i]][j,]<-unlist(lapply(1:length(String),function(x) {return(length(unlist(strsplit(SeqPFAM[[i]][j,x],split=""))))}))
     LengthSeqPFAM[[i]][j,which(LengthSeqPFAM[[i]][j,]==1)]<-0
    }
  
   if (j%%100==0) cat(j,", ",sep="")
   if (j%%1000==0) cat ("\n")
  }

# Compute the average sequence length of each PFAM along proteomes
MeanLengthSeqPFAM<-matrix(0,ncol=n.CorePFAM,nrow=n.proteomes)

for (j in 1:n.proteomes)
  {for (i in 1:n.CorePFAM)
    {Length<-LengthSeqPFAM[[i]][j,]
     Length<-Length[Length!=0]
     MeanLengthSeqPFAM[j,i]<-mean(Length)
    }
   
   if (j%%100==0) cat(j,", ",sep="")
   if (j%%1000==0) cat ("\n")
}


# Create txt file for each PFAM
PFAMi<-1 # Insert PFAM number
max.occ<-max(matrix_CorePFAM[PFAMi,])
listSeq<-t(matrix(unlist(SeqPFAM[PFAMi]), nrow=max.occ, byrow=T))
listProteins<-t(matrix(unlist(CorePROTEINS[PFAMi]), nrow=max.occ, byrow=T))
File<-c()
for(j in 1:n.proteomes)
  {NonZero<-which(listSeq[j,]!="0")
   IDProt.j<-paste(">",list_proteomes[j],"/",listProteins[j,], sep="")
   SeqProt.j<-paste(listSeq[j,],".", sep="")
   File<-c(File,rbind(IDProt.j[NonZero],SeqProt.j[NonZero]))
   
   if (j%%100==0) cat(j,", ",sep="")
   if (j%%1000==0) cat ("\n")
   
}

write.table(File, file = "SequencesCorePFAM_PFAMName.txt", quote=F,sep = "\t",
            row.names = FALSE, col.names = FALSE)

HeaderRow<-seq(1,nrow(File),2)
CoreSeqMat<-matrix(0,nrow=length(HeaderRow), ncol=3)
colnames(CoreSeqMat)<-c("Proteome","Protein_ID","Sequence")

for (i in 1:length(HeaderRow))
  {header<-unlist(strsplit(File[[1]][HeaderRows[i]], split="/"))
   CoreSeqMat[i,1]<-substr(header[1], start=2, stop=nchar(header)-4)
   CoreSeqMat[i,2]<-header[2]
   CoreSeqMat[i,3]<-substr(File[[1]][2*i], start=1, stop=nchar(File[[1]][2*i])-1)
  }

write.csv(CoreSeqMat, file="SequencesCorePFAM_PFAMName.csv", row.names=F, quote=F)
