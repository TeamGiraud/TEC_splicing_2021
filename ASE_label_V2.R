#### préparation du GTF
dat<-read.table(file="igenome_hg38.gtf",header=FALSE,sep="\t",colClasses=c('character','character','character','integer','integer','character','character','character','character'))

datex<-dat[dat[,3]=="exon",]
taille<-nrow(datex)



### extraction transcript_id
first<-regexpr("transcript_id",datex[1:taille,9])[1:taille]
ligne<-substr(datex[1:taille,9],first,nchar(datex[1:taille,9]))

fin<-regexpr(";",ligne)[1:taille]

transcript_id<-substr(ligne,nchar("transcript_id")+2,fin-1)



### extraction gene_id
first<-regexpr("gene_id",datex[1:taille,9])[1:taille]
ligne<-substr(datex[1:taille,9],first,nchar(datex[1:taille,9]))

fin<-regexpr(";",ligne)[1:taille]

gene_id<-substr(ligne,nchar("gene_id")+2,fin-1)



### addition de 2 variables pour comptage
gene_tag<-vector(length=taille,mode="integer")
gene_tag<-1
transcript_tag<-vector(length=taille,mode="integer")

datex2<-data.frame(datex,transcript_id,gene_id,gene_tag,transcript_tag)




### fichiers AS de MATS
#1
datSE<-read.table(file="fromGTF.SE.txt",header=TRUE,sep="\t",colClasses=c('integer','character','character','character','character','integer','integer','integer','integer','integer','integer'))

#2
datA3SS<-read.table(file="fromGTF.A3SS.txt",header=TRUE,sep="\t",colClasses=c('integer','character','character','character','character','integer','integer','integer','integer','integer','integer'))

#3
datA5SS<-read.table(file="fromGTF.A5SS.txt",header=TRUE,sep="\t",colClasses=c('integer','character','character','character','character','integer','integer','integer','integer','integer','integer'))

#4
datRI<-read.table(file="fromGTF.RI.txt",header=TRUE,sep="\t",colClasses=c('integer','character','character','character','character','integer','integer','integer','integer','integer','integer'))



### eviction trop petites différences

LIM<-12

datSE[,12]<-datSE[,7]-(datSE[,6]+1)
datSE<-datSE[datSE[,12]>=LIM,c(1:11)]

datA3SS[datA3SS[,"strand"]=="+",12]<-datA3SS[datA3SS[,"strand"]=="+",8]-(datA3SS[datA3SS[,"strand"]=="+",6])
datA3SS[datA3SS[,"strand"]=="-",12]<-datA3SS[datA3SS[,"strand"]=="-",7]-(datA3SS[datA3SS[,"strand"]=="-",9])
datA3SS<-datA3SS[datA3SS[,12]>=LIM,c(1:11)]

datA5SS[datA5SS[,"strand"]=="+",12]<-datA5SS[datA5SS[,"strand"]=="+",7]-(datA5SS[datA5SS[,"strand"]=="+",9])
datA5SS[datA5SS[,"strand"]=="-",12]<-datA5SS[datA5SS[,"strand"]=="-",8]-(datA5SS[datA5SS[,"strand"]=="-",6])
datA5SS<-datA5SS[datA5SS[,12]>=LIM,c(1:11)]

datRI[,12]<-datRI[,10]-(datRI[,9]+1)
datRI<-datRI[datRI[,12]>=LIM,c(1:11)]


### epurage "smooth": unique feature

#datSE[,12]<-datSE[,6]+datSE[,7]
#datSE<-datSE[match(unique(datSE[,12]),datSE[,12]),]

#datA3SS[,12]<-datA3SS[,6]+datA3SS[,7]
#datA3SS<-datA3SS[match(unique(datA3SS[,12]),datA3SS[,12]),]

#datA5SS[,12]<-datA5SS[,6]+datA5SS[,7]
#datA5SS<-datA5SS[match(unique(datA5SS[,12]),datA5SS[,12]),]

#datRI[,12]<-datRI[,6]+datRI[,7]
#datRI<-datRI[match(unique(datRI[,12]),datRI[,12]),]



### boucle sur datSE

datSE_taille<-nrow(datSE)
for (i in 1:datSE_taille){
	
	print(c(i,datSE_taille));
	
	datex2sub<-datex2[datex2[,"gene_id"]==datSE[i,2],]

	start0<-datSE[i,6]+1
	end0<-datSE[i,7]
	start_match0<-which(datex2sub[,4] %in% start0)
	end_match0<-which(datex2sub[,5] %in% end0)
	sel0<-as.vector(na.omit(start_match0[match(end_match0,start_match0)]))


	start1<-datSE[i,8]+1
	end1<-datSE[i,9]
	start_match1<-which(datex2sub[,4] %in% start1)
	end_match1<-which(datex2sub[,5] %in% end1)
	sel1<-as.vector(na.omit(start_match1[match(end_match1,start_match1)]))
	
	start2<-datSE[i,10]+1
	end2<-datSE[i,11]
	start_match2<-which(datex2sub[,4] %in% start2)
	end_match2<-which(datex2sub[,5] %in% end2)
	sel2<-as.vector(na.omit(start_match2[match(end_match2,start_match2)]))

	sel0NM<-as.vector(datex2sub[sel0,"transcript_id"])
	sel1NM<-as.vector(datex2sub[sel1,"transcript_id"])
	sel2NM<-as.vector(datex2sub[sel2,"transcript_id"])

	communNM_trans<-as.vector(na.omit(sel1NM[match(sel2NM,sel1NM)]))

	a<-as.vector(unique(datex2sub[datex2sub[,4]>end1 & datex2sub[,4]<start0,"transcript_id"]))
	b<-as.vector(unique(datex2sub[datex2sub[,5]>end1 & datex2sub[,5]<start0,"transcript_id"]))
	selEXCLU1<-as.vector(na.omit(b[match(a,b)]))

	a<-as.vector(unique(datex2sub[datex2sub[,4]>end0 & datex2sub[,4]<start2,"transcript_id"]))
	b<-as.vector(unique(datex2sub[datex2sub[,5]>end0 & datex2sub[,5]<start2,"transcript_id"]))
	selEXCLU2<-as.vector(na.omit(b[match(a,b)]))
	
	selEXCLU<-c(selEXCLU1,selEXCLU2)

	communNM<-communNM_trans[(communNM_trans %in% selEXCLU)==FALSE]


	incluNM<-as.vector(na.omit(sel0NM[match(communNM,sel0NM)]))
	excluNM<-communNM[(communNM %in% incluNM)==FALSE]

## ici astuce: match pour n'en prendre qu'un!!! pas mal...
	datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(incluNM,datex2[,"transcript_id"]),"gene_tag"]
	datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(excluNM,datex2[,"transcript_id"]),"gene_tag"]*9


	datex2[datex2[,"gene_id"]==datSE[i,2],"gene_tag"]<-as.vector(datex2[datex2[,"gene_id"]==datSE[i,2],"gene_tag"])[1]*10
}




### boucle sur datA3SS

datA3SS_taille<-nrow(datA3SS)
for (i in 1:datA3SS_taille){
	
	print(c(i,datA3SS_taille));
	
	datex2sub<-datex2[datex2[,"gene_id"]==datA3SS[i,2],]

	start0<-datA3SS[i,6]+1
	end0<-datA3SS[i,7]
	start_match0<-which(datex2sub[,4] %in% start0)
	end_match0<-which(datex2sub[,5] %in% end0)
	sel0<-as.vector(na.omit(start_match0[match(end_match0,start_match0)]))


	start1<-datA3SS[i,8]+1
	end1<-datA3SS[i,9]
	start_match1<-which(datex2sub[,4] %in% start1)
	end_match1<-which(datex2sub[,5] %in% end1)
	sel1<-as.vector(na.omit(start_match1[match(end_match1,start_match1)]))
	
	start2<-datA3SS[i,10]+1
	end2<-datA3SS[i,11]
	start_match2<-which(datex2sub[,4] %in% start2)
	end_match2<-which(datex2sub[,5] %in% end2)
	sel2<-as.vector(na.omit(start_match2[match(end_match2,start_match2)]))

	sel0NM<-as.vector(datex2sub[sel0,"transcript_id"])
	sel1NM<-as.vector(datex2sub[sel1,"transcript_id"])
	sel2NM<-as.vector(datex2sub[sel2,"transcript_id"])

	if(datex2sub[1,7]=="+"){
		a<-as.vector(unique(datex2sub[datex2sub[,4]>end2 & datex2sub[,4]<start0,"transcript_id"]))
		b<-as.vector(unique(datex2sub[datex2sub[,5]>end2 & datex2sub[,5]<start0,"transcript_id"]))
		selEXCLU<-as.vector(na.omit(b[match(a,b)]))
		}
	if(datex2sub[1,7]=="-"){
		a<-as.vector(unique(datex2sub[datex2sub[,4]>end0 & datex2sub[,4]<start2,"transcript_id"]))
		b<-as.vector(unique(datex2sub[datex2sub[,5]>end0 & datex2sub[,5]<start2,"transcript_id"]))
		selEXCLU<-as.vector(na.omit(b[match(a,b)]))
		}
	communNM<-sel2NM[(sel2NM %in% selEXCLU)==FALSE]
	
	incluNM<-as.vector(na.omit(sel0NM[match(communNM,sel0NM)]))
	excluNM<-as.vector(na.omit(sel1NM[match(communNM,sel1NM)]))

## ici astuce: match pour n'en prendre qu'un!!! pas mal...
	datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(incluNM,datex2[,"transcript_id"]),"gene_tag"]*2
	datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(excluNM,datex2[,"transcript_id"]),"gene_tag"]*9


	datex2[datex2[,"gene_id"]==datA3SS[i,2],"gene_tag"]<-as.vector(datex2[datex2[,"gene_id"]==datA3SS[i,2],"gene_tag"])[1]*10
}




### boucle sur datA5SS

datA5SS_taille<-nrow(datA5SS)
for (i in 1:datA5SS_taille){
	
	print(c(i,datA5SS_taille));
	
	datex2sub<-datex2[datex2[,"gene_id"]==datA5SS[i,2],]

	start0<-datA5SS[i,6]+1
	end0<-datA5SS[i,7]
	start_match0<-which(datex2sub[,4] %in% start0)
	end_match0<-which(datex2sub[,5] %in% end0)
	sel0<-as.vector(na.omit(start_match0[match(end_match0,start_match0)]))


	start1<-datA5SS[i,8]+1
	end1<-datA5SS[i,9]
	start_match1<-which(datex2sub[,4] %in% start1)
	end_match1<-which(datex2sub[,5] %in% end1)
	sel1<-as.vector(na.omit(start_match1[match(end_match1,start_match1)]))
	
	start2<-datA5SS[i,10]+1
	end2<-datA5SS[i,11]
	start_match2<-which(datex2sub[,4] %in% start2)
	end_match2<-which(datex2sub[,5] %in% end2)
	sel2<-as.vector(na.omit(start_match2[match(end_match2,start_match2)]))

	sel0NM<-as.vector(datex2sub[sel0,"transcript_id"])
	sel1NM<-as.vector(datex2sub[sel1,"transcript_id"])
	sel2NM<-as.vector(datex2sub[sel2,"transcript_id"])

	if(datex2sub[1,7]=="+"){
		a<-as.vector(unique(datex2sub[datex2sub[,4]>end0 & datex2sub[,4]<start2,"transcript_id"]))
		b<-as.vector(unique(datex2sub[datex2sub[,5]>end0 & datex2sub[,5]<start2,"transcript_id"]))
		selEXCLU<-as.vector(na.omit(b[match(a,b)]))
		}
	if(datex2sub[1,7]=="-"){
		a<-as.vector(unique(datex2sub[datex2sub[,4]>end2 & datex2sub[,4]<start0,"transcript_id"]))
		b<-as.vector(unique(datex2sub[datex2sub[,5]>end2 & datex2sub[,5]<start0,"transcript_id"]))
		selEXCLU<-as.vector(na.omit(b[match(a,b)]))
		}
	communNM<-sel2NM[(sel2NM %in% selEXCLU)==FALSE]

	incluNM<-as.vector(na.omit(sel0NM[match(communNM,sel0NM)]))
	excluNM<-as.vector(na.omit(sel1NM[match(communNM,sel1NM)]))

## ici astuce: match pour n'en prendre qu'un!!! pas mal...
	datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(incluNM,datex2[,"transcript_id"]),"gene_tag"]*3
	datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(excluNM,datex2[,"transcript_id"]),"gene_tag"]*9


	datex2[datex2[,"gene_id"]==datA5SS[i,2],"gene_tag"]<-as.vector(datex2[datex2[,"gene_id"]==datA5SS[i,2],"gene_tag"])[1]*10
}




### boucle sur datRI

datRI_taille<-nrow(datRI)
for (i in 1:datRI_taille){
	
	print(c(i,datRI_taille));
	
	datex2sub<-datex2[datex2[,"gene_id"]==datRI[i,2],]

	start0<-datRI[i,6]+1
	end0<-datRI[i,7]
	start_match0<-which(datex2sub[,4] %in% start0)
	end_match0<-which(datex2sub[,5] %in% end0)
	sel0<-as.vector(na.omit(start_match0[match(end_match0,start_match0)]))


	start1<-datRI[i,8]+1
	end1<-datRI[i,9]
	start_match1<-which(datex2sub[,4] %in% start1)
	end_match1<-which(datex2sub[,5] %in% end1)
	sel1<-as.vector(na.omit(start_match1[match(end_match1,start_match1)]))
	
	start2<-datRI[i,10]+1
	end2<-datRI[i,11]
	start_match2<-which(datex2sub[,4] %in% start2)
	end_match2<-which(datex2sub[,5] %in% end2)
	sel2<-as.vector(na.omit(start_match2[match(end_match2,start_match2)]))

	sel0NM<-as.vector(datex2sub[sel0,"transcript_id"])
	sel1NM<-as.vector(datex2sub[sel1,"transcript_id"])
	sel2NM<-as.vector(datex2sub[sel2,"transcript_id"])


	incluNM<-sel0NM
	excluNM<-as.vector(na.omit(sel1NM[match(sel2NM,sel1NM)]))

## ici astuce: match pour n'en prendre qu'un!!! pas mal...
	datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(incluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(incluNM,datex2[,"transcript_id"]),"gene_tag"]*4
	datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]<-datex2[match(excluNM,datex2[,"transcript_id"]),"transcript_tag"]+ datex2[match(excluNM,datex2[,"transcript_id"]),"gene_tag"]*9


	datex2[datex2[,"gene_id"]==datRI[i,2],"gene_tag"]<-as.vector(datex2[datex2[,"gene_id"]==datRI[i,2],"gene_tag"])[1]*10
}




### post-traitement

unique_NM<-unique(datex2[,"transcript_id"])
taille_unique<-length(unique_NM)

code<-vector(length<-taille_unique,mode="integer")
gene<-vector(length<-taille_unique,mode="character")

for (j in 1:taille_unique){
	gene[j]<-as.vector(datex2[datex2[,"transcript_id"]==as.vector(unique_NM[j]),"gene_id"])[1]
	code[j]<-sum(datex2[datex2[,"transcript_id"]==as.vector(unique_NM[j]),"transcript_tag"])
	print(c(j,taille_unique));
}



datfin<-data.frame(unique_NM,gene,code)

datfin2<-datfin[order(datfin[,"gene"],decreasing=F),]

write.csv(datfin2,file="ASE_label_V2_output.csv",row.names=FALSE)
