library(plyr)
library(ggplot2)
library(data.table)
#rm(list=ls())
pdfFile = 'effectPlots'
interest = c("x11","x90","x102")
#interest = TRUE

pooled = function(x,alpha=0.05,method='random'){
	if(method == 'random'){ #random effect
		Q = sum(x$w*x$g^2) - (sum(x$wf)^2)/sum(x$w)
		C = sum(x$w) - sum(x$w^2)/sum(x$w)
		tao2 = max((Q-nrow(x)+1)/C,0)
		v.r = x$se.g^2 + tao2
		w.r = 1/v.r
		pooled_g.r = sum(x$g*w.r)/sum(w.r) 
		pooled_se.r = sqrt(1/sum(w.r))
		pooled_CI.r = pooled_g.r + c(-1,1)*qnorm(1-alpha/2)*pooled_se.r
		ave_CI = mean(x$g) + c(-1,1)*qnorm(1-alpha/2)*sd(x$g)/sqrt(nrow(x))
		pval.r = 2*(1-pnorm(abs(pooled_g.r/pooled_se.r)))
		return(c(tao2=tao2,pooled_g=pooled_g.r,pooled_se=pooled_se.r,pval.z.r = pval.r,pooled_CI_lower=pooled_CI.r[1],pooled_CI_upper=pooled_CI.r[2],ave_CI_lower=ave_CI[1],ave_CI_upper=ave_CI[2]))
	} else{ #fixed effect
		pooled_g = sum(x$wf)/sum(x$w)
		pooled_se = sqrt(1/sum(x$w))
		pooled_CI = pooled_g + c(-1,1)*qnorm(1-alpha/2)*pooled_se
		ave_CI = mean(x$g) + c(-1,1)*qnorm(1-alpha/2)*sd(x$g)/sqrt(nrow(x))
	return(c(pooled_g=pooled_g,pooled_se=pooled_se,pooled_CI_lower=pooled_CI[1],pooled_CI_upper=pooled_CI[2],ave_CI_lower=ave_CI[1],ave_CI_upper=ave_CI[2]))
		
		}
}

fisher = function(pval){
	chisq = -2*sum(log(pval))
	p.meta = 1-pchisq(chisq,df = 2*length(pval))
	return(c(chisq.stat = chisq,p.meta = p.meta))
}

importFile = function(path,pattern,subFileName){
	filenames = list.files(path=path,pattern = pattern)
	listData = llply(filenames,read.delim,header=T)
	filenames = gsub(subFileName,'',filenames)
	names(listData) = filenames
	return(rbindlist(listData))	
}

pvalFile = importFile(path = '.', pattern= '.txt', subFileName = '.txt') ##import stat test file

rownames(pvalFile) = make.names(pvalFile$Dataset,unique=T)
pvalFile = pvalFile[,-1]
pvalFile = pvalFile[pvalFile$name%in%interest,]
if(isTRUE(interest)){
	queries = unique(pvalFile$name)
} else{
	queries = interest
}
rm(interest)

effectPlot = function(query, effectSizes){
	es_query = effectSizes[effectSizes$query.Symbol == query,] ## effect size info for query
	es_pooled = pooled(x = es_query)
	meta_chisq = fisher(es_query$FDRt2)
	return(list(es_query, es_pooled, meta_chisq))
}

plotES = function(pooled, es, query = query){
	esPlot = ggplot() + geom_polygon(data=pooled,aes(x=x,y=y),fill="#FFCC33" + geom_point(data=es,aes(x=g,y=code,size=w),pch=15) + geom_vline(xintercept=0,lty=2) + theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_rect(colour='black',size=1),legend.position='none') + ylab("") + xlab("Standardized Mean Difference (log2 scale)") + ggtitle(sprintf('%s',query)) + geom_errorbarh(data=a,aes(x=g,xmin=l,xmax=u,y=code),height=0,alpha=0.5,size=0.7)+scale_y_continuous(breaks=c(a$code,0,-1),labels=c(samp,'Summary','Summary_average'))
	return(esPlot)
}

pdf(sprintf('%s.pdf',pdfFile))
for(k in 1:length(queries)){
	query = queries[k]
	use = effectPlot(effectSizes = pvalFile)
	fooSet = use[[1]] ## effect size info for the foo query
	pooled_ES = use[[2]]
	meta.chisq = use[[3]]
	study = data.frame(g=fooSet$g,l=fooSet$CI.lower,u=fooSet$CI.upper,code=seq(1,nrow(fooSet)),rownames(fooSet),w=1/(1/fooSet$w+pooled_ES[1]))

	pooled_g = use[[2]]$pooled_g
	
	pooled = data.frame(y=c(0,0.2,0,-0.2),x=c(use[[2]]$'pooled_CI_lower',pooled_g,use[[2]]$'pooled_CI_upper',pooled_g))
	plotES(pooled = pooled, es = study)
	write.table(t(c(query,pooled_ES,meta.chisq)),'pooledES',append=(i>1),col.names=(i == 1),quote=F,sep='\t')
}
dev.off()
