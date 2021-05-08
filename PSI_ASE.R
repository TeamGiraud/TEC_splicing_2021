#### load gene frequency file from Tophat2
# here example: genes_mTEChi.fpkm_tracking
dat<-read.table(file="genes_mTEChi.fpkm_tracking", header=TRUE,sep="\t")

#### load isoform frequency file from Tophat2
# here example: isoforms_mTEChi.fpkm_tracking
dat_iso<-read.table(file="isoforms_mTEChi.fpkm_tracking", header=TRUE,sep="\t")

#### load ASE_label_output.csv file from ASE_label_V2.R or ASE_label_V3.R
# here example: ASE_label_output.csv
dat_incl<-read.csv(file="ASE_label_output.csv",header=TRUE)

dat<-dat[dat[,"gene_short_name"]!="-",]

dat<-dat[dat$FPKM>0,]






I<-vector(mode="numeric",length=nrow(dat))

glob<-matrix(data=NA,nrow=1,ncol=2)
colnames(glob)<-c("X1","X2")

for (i in 1:nrow(dat)){
	sub_incl<-dat_incl[dat_incl$gene==as.character(dat[i,"gene_id"]),]
	if (nrow(sub_incl)==0|nrow(sub_incl)==1){I[i]<-NA} else {
		MAX<-nchar(as.character(max(sub_incl$code)))
	
		Isub<-vector(mode="numeric",length=MAX)
		
		for (k in 1:MAX){
		

			comp<-as.integer(substr(as.character(sub_incl$code),nchar(as.character(sub_incl$code))-k+1,nchar(as.character(sub_incl$code))-k+1))
			comp[is.na(comp)]<-0
			
			if (sum(comp)/length(comp[comp!=0])!=9 & sort(comp==9,decreasing=T)[1]==TRUE){
				NM_incl<-as.vector(sub_incl[comp!=0&comp!=9,"unique_NM"])
				NM_comp<-as.vector(sub_incl[comp!=0,"unique_NM"])
				
				value1_incl<-sum(dat_iso[match(NM_incl,dat_iso$tracking_id),"FPKM"])
				value1_comp<-sum(dat_iso[match(NM_comp,dat_iso$tracking_id),"FPKM"])
				Isub[k]<-value1_incl/value1_comp
				temp<-data.frame(dat[i,"gene_id"],Isub[k])
				colnames(temp)<-c("X1","X2")
			}else 
				{Isub[k]<-NA
				temp<-c(NA,NA)}

			glob<-rbind(data.frame(glob,row.names=NULL),temp)

		}
		

		
		I[i]<-mean(Isub,rm.na=TRUE)
		
	print(c(i,nrow(dat)))

}
}

dat2<-data.frame(dat,I)


write.csv(na.omit(glob),file="PSI_ASE_output.csv",row.names=FALSE)
