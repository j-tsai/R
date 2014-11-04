#!/usr/bin/Rscript
#eg.Rscript ./DESeq.R ./ rawCounts.txt design_trial.txt effect1 testing local
args = commandArgs(TRUE)

#Anders S and Huber W (2010). “Differential expression analysis for sequence count data.” Genome Biology, 11, pp. R106. http://dx.doi.org/10.1186/gb-2010-11-10-r106, http://genomebiology.com/2010/11/10/R106/
if(length(args)<1)
	args = c("--help")
if("--help" %in% args){
	cat("
		Arguments:
		dir : output directiory
		raw_counts : input raw counts file
		design : design - samples in rows and conditions in columns
		effect : variable in the design matrix of interest
		o : output file name
		fitType (optional) : fit type for estimating dispersion
		
		--help
		output: o_effect.txt p-values, adjusted p-values, normalized counts
		Example:
		Rscript ./DESeq_1014.R ./ rawCounts_eg.txt design_eg.txt effect1 testing local \n\n")
	q(save='no')
}

if("DESeq" %in% rownames(installed.packages()) == FALSE) {
	source("http://bioconductor.org/biocLite.R")
	biocLite("DESeq")}
library("DESeq")

write.t = function(data,outfile){
	dt = data.frame(GeneSymbol = rownames(data),data)
	write.table(dt,sprintf('%s.txt',outfile),sep='\t',quote=F,row.names=F)
}


setwd(args[1])

raw = read.delim(args[2],row.names=1)
rawDesign = read.delim(args[3],row.names=1)

effect = args[4]


num_var = ncol(rawDesign)

if(num_var==1){
	rawDesign=rawDesign[,1]
}else if(num_var>2){
	print('#variables > 2. Abort')
	q(save='no')
}
	

	

cds = newCountDataSet(raw, rawDesign) 

cds = estimateSizeFactors(cds)

print('size factors:')
sizeFactors(cds)


if(is.na(args[6])) 
	args[6]='parametric'

	
cds = estimateDispersions(cds,fitType=args[6])


png(sprintf('dispPlot_%s.png',args[5]))

plotDispEsts( cds )
dev.off()

write.t(fData(cds),sprintf('dispersions_%s',args[5]))
write.t(counts(cds,normalized=TRUE),sprintf('normalizedCounts_%s',args[5]))
write.table(sizeFactors(cds),sprintf('sizeFactors_%s',args[5]))

pData(cds)
ncol(pData(cds))


ptm <- proc.time()
if(num_var==1){
	fit1 = fitNbinomGLMs( cds, count ~ rawDesign)
	fit0 = fitNbinomGLMs(cds, count~ 1)
} else{
	fit1 = fitNbinomGLMs( cds, count ~ rawDesign[,colnames(rawDesign)!=args[4]]+rawDesign[,args[4]])
	fit0 = fitNbinomGLMs(cds, count~ rawDesign[,colnames(rawDesign)!=args[4]])
} 
proc.time() - ptm


sprintf('#converged = %s',sum(fit0$converged ==TRUE &fit1$converged ==TRUE ))


pvalsGLM = nbinomGLMTest(fit1,fit0) 

names(pvalsGLM) = featureNames(cds)


padj = p.adjust(pvalsGLM,method='BH')

sprintf("# genes adjusted p-value < 0.05: %s",sum(padj< 0.05, na.rm=TRUE))


write.t(data.frame(rawPval = pvalsGLM,fdr=padj,counts(cds,normalized=T)),sprintf('%s_%s',args[5],args[4]))