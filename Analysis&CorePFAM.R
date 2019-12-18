matrix_complete<-read.csv(file="matrix_PFAMvsPROTEOMES_bacteria.csv",header=F,sep=";")
matrix_ab<-as.matrix(matrix_complete[-c(1),-c(1)])
rownames(matrix_ab) <-NULL
colnames(matrix_ab) <-NULL
class(matrix_ab) <- "numeric" 
rm(matrix_complete)

n.proteomes<-ncol(matrix_ab)
n.PFAM<-nrow(matrix_ab)

af<-rowSums(matrix_ab)
lp<-colSums(matrix_ab)

matrix.presence<-matrix_ab
matrix.presence[matrix.presence>0]<-1

presence.PFAMinProteomes<-rowSums(matrix.presence)

MaxOccInProteomes<-rep(0,n.PFAM)
for (i in 1:n.PFAM) MaxOccInProteomes[i]<-max(matrix_ab[i,1:n.proteomes])

MaxOcc<-4 # Set maximal occurrence to define a Core PFAM
rarePFAM<-which(MaxOccInProteomes<=MaxOcc)
presence.rarePFAMinProteomes<-presence.PFAMinProteomes[rarePFAM]

percentagePROT<-seq(0,100,by=1)
rarePFAMvsPercentage<-rep(0,length(percentagePROT))
for (i in 1:length(percentagePROT))
rarePFAMvsPercentage[i]<-length(which(presence.rarePFAMinProteomes>((n.proteomes)*percentagePROT[i]/100)))

#################  FIND CORE-PFAMs #################
# 1: Set the wished percentage of minimal prevalence
MinPrev<-0.95

# 2: Select the corresponding PFAM
plus<-which(presence.PFAMinProteomes>round(n.proteomes*MinPrev))

# 3: The core PFAM can be obtained by intersecting the most prevalent and the rarest PFAM
sel<-unique(intersect(plus,rarePFAM))

# 4: Check whether one or more proteomes have been missed
sel_mat<-matrix_ab[sel,]
sel_vec<-colSums(sel_mat)

missed_prot<-which(sel_vec==0)

# If no proteomes have been missed, go to step 8. Otherwise proced to step 5:

# 5: Check whether there is a PFAM covering all missed proteomes. If so go to 6. Otherwise go to step 7.
length(intersect(rarePFAM,which(rowSums(missed_prot_mat)==length(which(sel_vec==0)))))

#6: Add the PFAM to the list of CorePFAM
z<-intersect(rarePFAM,which(rowSums(missed_prot_mat)==length(which(sel_vec==0))))
max.z<-max(presence.PFAMinProteomes[z]/n.proteomes*100)
index.PFAM.to.add<-which(presence.PFAMinProteomes[z]/n.proteomes*100==max.z)
sel<-c(sel,z[index.PFAM.to.add])

# 7: Find the most prevalent PFAM needed to cover the missed proteomes
PFAM.to.add<-rep(0,length(missed_prot))

for (i in 1:length(missed_prot))
 {missed_prot_i<-matrix.presence[,missed_prot[i]]
  z<-intersect(rarePFAM,which(missed_prot_i==1))
  max.z<-max(presence.PFAMinProteomes[z]/n.proteomes*100)
  index.PFAM.to.add<-which(presence.PFAMinProteomes[z]/n.proteomes*100==max.z)
  PFAM.to.add[i]<-z[index.PFAM.to.add]}
sel<-c(sel,PFAM.to.add)

# 8: Reduce Core PFAM number, if possible
CorePFAM<-c(sel[c(22,23,31,33,34,36,38,39)],PFAM.to.add)

# 9: Save the list of Core PFAM
save(CorePFAM,file="list_CorePFAM_bacteria")
