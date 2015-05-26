library("ggplot2")
write.table0 = function(data,outfile){
	write.table(data,file = sprintf('%s.txt',outfile),quote=F,sep='\t',row.names=F)
}

		# dir : output directiory
		# raw_counts : input raw counts file
		# ctrl = control
		# t = target
		# output : output file name

dir = "~/Dropbox/qPCR_script/output3"
output = 'outputfile' 

raw_counts = '~/Dropbox/qPCR_script/efsun_qpcr.txt'
ctrl = "EndoC-SIX3"
t = 'ACTB'

sample0 = NULL
target0 = NULL

setwd(dir)
x = read.delim(raw_counts,header=T)
x = data.frame(name = paste(x$Sample.Name,x$Target.Name,sep='::'),x)
sample = unique(x$S); print(sprintf('Samples: %s', paste(sample,collapse=' | ')))

targets = as.character(unique(x$Target.Name[x$Sample.Name==ctrl])); print(sprintf('Targets in ref: %s', paste(targets,collapse=' | ')))

s0 = function(x){
	return(c(mean(x,na.rm=T),sd(x,na.rm=T))) 
}


ag = aggregate(as.numeric(as.character(x$C)),by = list(x$name),s0)
tmp = t(as.data.frame(strsplit(as.character(ag$Group.1),'::')))
ag = data.frame(sample = tmp[,1],target = tmp[,2],gr = ag$Group.1,mean = ag$x[,1],sd = ag$x[,2])
remove(tmp)
### dCt
s1 = NULL; s2 = NULL; s3 = NULL
for(i in sample){
	ref = ag$mean[ag$sample == i & ag$target == t]
	for(j in unique(ag$target)){
		if(length(ag$mean[ag$sample == i & ag$target == j])>0){
			s1 = append(s1,i)
			s2 = append(s2,j)
			s3 = append(s3,ag$mean[ag$sample == i & ag$target == j] - ref)
		}
	}
}
dct = data.frame(name = paste(s1,s2,sep="::"),sample = s1, target = s2, dCt = s3, FC2dCt = 2^-s3) ### one less SIX2
### ddCt
s1 = NULL; s2 = NULL; s3 = NULL

for(i in targets){
	ref = dct$dCt[dct$sample == ctrl & dct$target == i]
	for(j in sample){
		if(length(dct$dCt[dct$sample == j & dct$target == i])>0){
			s1 = append(s1,j)
			s2 = append(s2,i)
			s3 = append(s3,dct$dCt[dct$sample == j & dct$target == i] - ref)
		}
	}
}
ddct = data.frame(name = paste(s1,s2,sep="::"),sample = s1, target = s2, ddCt = s3, FC2ddct = 1/2^(s3)) 
dct = dct[,-c(2,3)]
ct = merge(dct,ddct,by='name',all=T); ct = ct[,-c(4,5)]
rm(ddct);rm(dct)
ct = merge(ag,ct,by.x='gr',by.y = 'name',all=T); ct = ct[,-c(2,3)]
rm(ag)
xx = merge(x,ct,by.x='name',by.y='gr')
rm(ct)
xx = data.frame(xx,Ct_mean35 = mapply(function(x,y) if(!is.na(x) & x > 35){
	y
} else{
	''
},as.numeric(as.character(xx$C)),as.character(xx$Well)),
Ct_sd1 = sapply(xx$sd,function(x) if(!is.na(x) & x > 1){
	'**'
} else{
	''
})
)
write.table0(xx,sprintf('%s_aggregate',output))

xx = merge(aggregate(xx[c(3,4,6,7,8,9,10,11,13)],xx[1],unique),aggregate(xx[12],xx[1],paste,collapse=' '),by.x='name',by.y='name')
write.table0(xx,sprintf('%s_collapse',output))

### plot dCt
if(is.null(sample0))
	sample0 = sample
if(is.null(target0))
	target0 = targets[targets!=t]

for(i in sample0){
	c = ggplot(xx[xx$Sample.Name==i,],aes(x=Target.Name,y=mean,width=.5))+theme_bw()+theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
	#limits = aes(ymax = mean + sd, ymin = mean - sd)
	png(sprintf('%s_%s_dCt.png',output,i))
	print(c + geom_bar(stat='identity',position='dodge')+xlab('Targets')+ylab('dCt')+ggtitle(i))
	#+geom_errorbar(limits,width=0.25)
	dev.off()
}

for(i in target0){
	if(sum(is.nan(xx[xx$Target.Name==i,'mean'])) != sum(xx$Target.Name==i)){
		c = ggplot(xx[xx$Target.Name==i,],aes(x=Sample.Name,y=FC2dCt,width=.5))+theme_bw()
		png(sprintf('%s_%s_FC2dCt.png',output,i))
		print(c + geom_bar(stat='identity',position='dodge')+xlab('Targets')+ylab('dCt')+ggtitle(i))
	dev.off()

	}
}
