##Rscript: completeTCR.R (chiotti
##Merge the 'long' ids produced by the modified MiTCR softeware (CompleteTCR) and stored in results.txt with those in the new descToID_mapping.txt file to get the input read names; match up read ids between alpha and beta inputs; get complete clone frequencies and percent reads
##INPUT: 
##OUTPUT: 

library(plyr)

args=commandArgs(TRUE)

outDir=args[1]
resultFileA=args[2]
mapFileA=args[3]
resultFileB=args[4]
mapFileB=args[5]
name=args[6]

outdir=gsub("/","//",outDir)
resultA=gsub("/","//",resultFileA)
mapfileA=gsub("/","//",mapFileA)
resultB=gsub("/","//",resultFileB)
mapfileB=gsub("/","//",mapFileB)

rsltA=read.table(file=resultA,header=T,row.names=NULL,sep="\t",skip=1)
mapA=read.table(file=mapfileA,header=F,row.names=NULL,sep="\t")
colnames(mapA)=c('LongID','readID','readSeq')
rsltB=read.table(file=resultB,header=T,row.names=NULL,sep="\t",skip=1)
mapB=read.table(file=mapfileB,header=F,row.names=NULL,sep="\t")
colnames(mapB)=c('LongID','readID','readSeq')

mappedA=paste(outdir,"//",name,".mappedAlpha.txt", sep="")
mappedB=paste(outdir,"//",name,".mappedBeta.txt",sep="")
pairedOut=paste(outdir,"//",name,".alphaBetaPairs.txt",sep="")
aggPairs=paste(outdir,"//",name,".aggregatedPairs.txt",sep="")

matchIDs=function(x){
#	wkg=data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("VAlleles", "DAllele","JAlleles","ntCDR3","aaCDR3"))))
	wkg=data.frame(matrix(vector(), 0, 12, dimnames=list(c(), c("VAlleles", "DAllele","JAlleles","ntCDR3","aaCDR3","lastV","firstD","lastD","firstJ","VDinserts","DJinserts","Total_inserts"))))
	for (n in 1:nrow(x)) {
		idstring=as.character(x[n,'IDs'])
		ids=as.data.frame(strsplit(idstring,split=";"))
		colnames(ids)="ID"
#		cnst=c(as.character(x[n,'V.alleles']),as.character(x[n,'D.alleles']),as.character(x[n,'J.alleles']),as.character(x[n,'CDR3.nucleotide.sequence']),as.character(x[n,'CDR3.amino.acid.sequence']))
		cnst=c(as.character(x[n,'V.alleles']),as.character(x[n,'D.alleles']),as.character(x[n,'J.alleles']),as.character(x[n,'CDR3.nucleotide.sequence']),as.character(x[n,'CDR3.amino.acid.sequence']),as.character(x[n,'Last.V.nucleotide.position']),as.character(x[n,'First.D.nucleotide.position']),as.character(x[n,'Last.D.nucleotide.position']),as.character(x[n,'First.J.nucleotide.position']),as.character(x[n,'VD.insertions']),as.character(x[n,'DJ.insertions']),as.character(x[n,'Total.insertions']))
		con=t(as.data.frame(rep(as.data.frame(cnst),nrow(ids))))
		join=cbind(ids,con)
		wkg=rbind(wkg,join)
	}
#	colnames(wkg)=c("LongID","VAlleles","DAllele","JAlleles","ntCDR3","aaCDR3")	
	colnames(wkg)=c("LongID","VAlleles","DAllele","JAlleles","ntCDR3","aaCDR3","lastV","firstD","lastD","firstJ","VDinserts","DJinserts","Total_inserts")
	return(wkg)
}

wkgA=matchIDs(rsltA)
wkgB=matchIDs(rsltB)

write.table(wkgA, file=mappedA,col.names=T,row.names=F,sep="\t",quote=F)
write.table(wkgB, file=mappedB,col.names=T,row.names=F,sep="\t",quote=F)

mrgA=merge(mapA[,1:2],wkgA,by='LongID',all=T)
mrgB=merge(mapB[,1:2],wkgB,by='LongID',all=T)
mrgIDs=merge(mrgA[,-1],mrgB[,-1],by='readID',all=T)

pairs=mrgIDs[,-3]
#colnames(pairs)=c('readID','alpha_VAlleles','alpha_JAlleles','alpha_ntCDR3','alpha_aaCDR3','beta_VAlleles','beta_DAllele','beta_JAlleles','beta_ntCDR3','beta_aaCDR3')
colnames(pairs)=c('readID','alpha_VAlleles','alpha_JAlleles','alpha_ntCDR3','alpha_aaCDR3','alpha_lastV','alpha_firstD','alpha_lastD','alpha_firstJ','alpha_VDinserts','alpha_DJinserts','alpha_Total_inserts','beta_VAlleles','beta_DAllele','beta_JAlleles','beta_ntCDR3','beta_aaCDR3','beta_lastV','beta_firstD','beta_lastD','beta_firstJ','beta_VDinserts','beta_DJinserts','beta_Total_inserts')

truePairs=na.omit(pairs)
#colnames(truePairs)=c('readID','alpha_VAlleles','alpha_JAlleles','alpha_ntCDR3','alpha_aaCDR3','beta_VAlleles','beta_DAllele','beta_JAlleles','beta_ntCDR3','beta_aaCDR3')
colnames(truePairs)=c('readID','alpha_VAlleles','alpha_JAlleles','alpha_ntCDR3','alpha_aaCDR3','alpha_lastV','alpha_firstD','alpha_lastD','alpha_firstJ','alpha_VDinserts','alpha_DJinserts','alpha_Total_inserts','beta_VAlleles','beta_DAllele','beta_JAlleles','beta_ntCDR3','beta_aaCDR3','beta_lastV','beta_firstD','beta_lastD','beta_firstJ','beta_VDinserts','beta_DJinserts','beta_Total_inserts')

#write.table(pairs,file=pairedOut,col.names=T,row.names=F,sep="\t",quote=F)
#write.table(truePairs,file=pairedOut,col.names=T,row.names=F,sep="\t",quote=F)

#name=""
#fx=function(name){
#outdir=""
#infile=paste(outdir,"//",name,"//",name,".alphaBetaPairs.txt",sep="")
#truePairs=read.table(file=infile,header=T,row.names=NULL,sep="\t")
#aggPairs=paste(outdir,"//",name,".aggregatedPairs.txt",sep="")

#agg=count(truePairs,c('alpha_VAlleles','alpha_JAlleles','alpha_ntCDR3','alpha_aaCDR3','beta_VAlleles','beta_DAllele','beta_JAlleles','beta_ntCDR3','beta_aaCDR3'))
agg=count(truePairs,c('alpha_VAlleles','alpha_JAlleles','alpha_ntCDR3','alpha_aaCDR3','alpha_lastV','alpha_firstD','alpha_lastD','alpha_firstJ','alpha_VDinserts','alpha_DJinserts','alpha_Total_inserts','beta_VAlleles','beta_DAllele','beta_JAlleles','beta_ntCDR3','beta_aaCDR3','beta_lastV','beta_firstD','beta_lastD','beta_firstJ','beta_VDinserts','beta_DJinserts','beta_Total_inserts'))
agg=agg[order(agg['freq'],decreasing=T),c(1:ncol(agg))]

pcnt=numeric(length=nrow(agg))

for (i in 1:nrow(agg)){
	pcnt[i]=round(agg[i,10]/nrow(truePairs),9)
}

pcnt=as.data.frame(pcnt)
agg=cbind(agg,pcnt)

write.table(agg,file=aggPairs,col.names=T,row.names=F,sep="\t",quote=F)
#}
q()
