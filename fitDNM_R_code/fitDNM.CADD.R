library(foreach)
library(iterators)
library(doParallel)

nodes =detectCores()
registerDoParallel(nodes)



args <-commandArgs(trailingOnly=TRUE)


source(args[7])

nmale = as.integer(args[1])
nfemale = as.integer(args[2])
DNMFile = args[3]
gene_mu_rate_file = args[4]
gene_polyphen_file = args[5]
resultFile = args[6]
annotation = toString(args[8])

DNM = read.table(DNMFile,header=T,stringsAsFactors=F,colClasses=c("character","integer","character","character","character"))
#DNM = read.table(DNMFile,header=T,stringsAsFactors=F)

if (nrow(DNM) > 0){
	gene_mu_rate = read.table(gene_mu_rate_file,header=T,stringsAsFactors=F)
	gene_polyphen = read.table(gene_polyphen_file,stringsAsFactors=F,header=T,fill=T)
	#gene_polyphen_new = gene_polyphen[-na.list,]
	#gene_mu_rate_new = gene_mu_rate[-na.list,]



	na.list = which(gene_polyphen[,5] < 0|gene_polyphen[,6] < 0|gene_polyphen[,7] < 0|gene_polyphen[,8] < 0|is.na(gene_polyphen[,5]) |is.na(gene_polyphen[,6])|is.na(gene_polyphen[,7])|is.na(gene_polyphen[,8]))

	gene_polyphen_new = gene_polyphen
	#gene_polyphen_new <- gene_polyphen_new[-1,]
	#gene_polyphen_new <- gene_polyphen_new[-nrow(gene_polyphen_new),]

	gene_mu_rate_new = gene_mu_rate
	#gene_polyphen_new = gene_polyphen[-na.list,]
	#gene_mu_rate_new = gene_mu_rate[-na.list,]

	geneIDs = as.character(unique(DNM[,5]))
	ngene = length(geneIDs)



	##convert allele to number A->1;T->2;C->3;G->4
	allele_2_num = function(allele){
		if(allele=="A"){
			return(1)
		}else if(allele=="T"){
			return(2)
		}else if(allele=="C"){
			return(3)
		}else if(allele=="G"){
			return(4)
		}
	}



	compute_pvalue = function(igene){
	  DNM_tmp = DNM[which(DNM[,5] ==geneIDs[igene]),]
	  if(DNM_tmp[1,1] =="X"){
	  	nsample = (nmale+nfemale*2)/2
	  }else{
	  	nsample = nmale+nfemale
	  }

	  nsnv_o =length(which(gene_polyphen[,1]==geneIDs[igene]))
	  gene_pos = gene_polyphen_new[which(gene_polyphen_new[,1]==geneIDs[igene]),3]
	  gene_polyphen_tmp_1 = as.matrix(gene_polyphen_new[which(gene_polyphen_new[,1]==geneIDs[igene]),5:8])
	  gene_mu_rate_tmp = as.matrix(gene_mu_rate_new[which(gene_mu_rate_new[,1]==geneIDs[igene]),5:8])
	  gene_reference = gene_polyphen_new[which(gene_polyphen_new[,1]==geneIDs[igene]),4]

	  ##recode for LOF mutations
	  gene_polyphen_tmp = gene_polyphen_tmp_1
	  #gene_polyphen_tmp[gene_polyphen_tmp_1==2] = 1
	  #gene_polyphen_tmp[gene_polyphen_tmp_1==3] = 0

	 nsnv = dim(gene_polyphen_tmp)[1]

	 stat_weights = gene_polyphen_tmp - matrix(rep(apply((gene_polyphen_tmp* gene_mu_rate_tmp),1,sum),4),ncol=4)
	  exp_lambda = nsample*2*gene_mu_rate_tmp

	  if(any(is.na(gene_polyphen_tmp))){
	  	cat(igene,"NA weights \n")
	  }
	    nsnv = dim(gene_polyphen_tmp)[1]




	  ndenovo = dim(DNM_tmp)[1]
	  scores=0
	  scores2=0 ##for unweighted

	   if(dim(gene_mu_rate_tmp)[1] <=0){
	  		pvalue =1
	  		pvalue.unweight =1
	  		return(c(geneIDs[igene],nsample,nsnv_o,nsnv,ndenovo,scores2,round(scores,3),pvalue,pvalue.unweight,NA))
	  }
	  	ref_pos = which(gene_mu_rate_tmp==0,arr.ind=T)
	  	ref_pos_v = apply(ref_pos,1,function(x) ((x[2]-1)*nsnv+x[1]))


	  if(ndenovo ==0 ) {
	  	return(c(geneIDs[igene],nsample,nsnv,ndenovo,1,NA,1,1,NA))
	  }
	  for(idenovo in c(1:ndenovo)){
	  	###check if the position exist in the ccds(not included promoter region)
	  	if(DNM_tmp[idenovo,2] %in% gene_pos){
	  		if(DNM_tmp[idenovo,3] != gene_reference[which(gene_pos==DNM_tmp[idenovo,2])]){
	  			cat(DNM_tmp[idenovo,1],DNM_tmp[idenovo,2],DNM_tmp[idenovo,3],"donot match ref\n")
	  		}
	        scores =scores + stat_weights[which(gene_pos==DNM_tmp[idenovo,2]), allele_2_num(DNM_tmp[idenovo,4])]
	  		scores2 = scores2 +1
	  	}
	  }
	    pvalue = double_saddle_point_approximation(scores,c(exp_lambda)[-ref_pos_v],c(stat_weights)[-ref_pos_v])

		##compute the pvalue with poisson(unweighted)
		pvalue.unweight = 1-ppois(scores2-1,sum(exp_lambda))



		cat(c(geneIDs[igene],nsample,nsnv_o,nsnv,ndenovo,round(scores,3),pvalue,pvalue.unweight,"\n"))
		return(c(geneIDs[igene],nsample,nsnv_o,nsnv,ndenovo,round(scores,3),pvalue,pvalue.unweight))
	}


		stat.list = foreach(igene = 1:ngene, .combine=rbind) %dopar%{

			compute_pvalue(igene)
		}

	if(typeof(stat.list) == "character"){stat.list <- as.data.frame(t(stat.list))}

		colnames(stat.list) = c("geneID","nsample.adj","nsnv","nsnv.analysis","ndenovo","score","pvalue.fitDNM","pvalue.Poisson")

	    write.table(stat.list,file=resultFile,quote=F,row.names=F,col.names=T)
}else{
	out_file <- data.frame(c(paste("No SNVs in",annotation,sep =" ")))
	print(annotation)
	colnames(out_file) <- c("header")
	write.table(out_file,file=resultFile,quote=F,row.names=F,col.names=T)
}
