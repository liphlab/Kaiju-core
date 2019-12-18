list_proteomes<-read.table(file="list_proteomes_bacteria.txt",sep=" ",header=F,as.is=TRUE)
list_PFAM<-read.table(file="list_PFAM_bac.txt",sep=" ",header=F,as.is=TRUE)

n_prot<-ncol(list_proteomes)
n_PFAM<-nrow(list_PFAM)

matrix_ab<-matrix(0,nrow=n_PFAM,ncol=n_prot)

for (i in 1:n_prot)
  {name.prot<-list_proteomes[i]
   prot_list<-read.table(file=paste("abundances_proteomes_bacteria/",name.prot,sep=""),header=FALSE,as.is=TRUE)
   PFAM_involved<-match(prot_list[,1],list_PFAM[,1])
   matrix_ab[PFAM_involved,i]<-prot_list[,2]
   if (i%%10==0) cat(i," ")
   if (i%%100==0)cat("\n")
  }

matrix_complete<-matrix(0,nrow=n_PFAM+1,ncol=n_prot+1)
matrix_complete[1,1]<-"PFAMvsPROTEOMES"
matrix_complete[1,2:(n_prot+1)]<-paste(list_proteomes)
matrix_complete[2:(n_PFAM+1),1]<-list_PFAM[,1]
matrix_complete[2:(n_PFAM+1),2:(n_prot+1)]<-matrix_ab

#write.table(matrix_complete,file="matrix_PFAMvsPROTEOMES_bacteria.csv",row.names=F,col.names=F,sep=";")
