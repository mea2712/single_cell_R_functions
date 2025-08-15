
getAnn_f <- function(dataset,version){
# dataset: biomart dataset
# version: ensembl version
  att = c("ensembl_gene_id","ensembl_gene_id_version","external_gene_name","ensembl_transcript_id_version") # Can add more, but these are the minimal
  # Get IDs
  mrt = useEnsembl(biomart='genes',dataset=dataset,version=version)
  gid = getBM(attributes=att, mart=mrt)
  # Remove ensembl_gene_id that match to multiple external_gene_names
  tmp = tibble(gid) %>% dplyr::group_by(ensembl_gene_id) %>% summarize(n=n_distinct(external_gene_name)) %>% dplyr::filter(n>1) %>% dplyr::select(ensembl_gene_id) %>% unlist()
  # Sometimes, some cleaning needed
  if(length(tmp)>0) {gid = gid %>% dplyr::filter(! ensembl_gene_id %in% tmp); rm(tmp) }
  idx= which(gid$ensembl_gene_id_version==gid$ensembl_transcript_id_version)
  if(length(idx)>0) {gid=gid[-idx,] ; rm(idx)}
  gid=unique(gid)
  return(gid)
}

procAnn_f <- function(args1, args2, args3, args4,args5){
	
# Get working annotation - needs server with internet connection -
# args5=as.character('format') - can take value 'download' (download from AnnHub), 'format' (format an existing gtf), 'nothing' (do nothing, just subset to require fields)
# args1=as.character('path_to_gtf') can take value 'path_to_gtf' (if args0  is 'format' or 'nothing') 'NULL' if args0 is download
# args2=c(see below)
# args3=as.character(no) - or 'yes', whether you'd like to write and rds file with the annotation to disk
# args4=as.character('path_to_output_file') - 'NULL' if args3 is no

# gtf must have fields for ensembl_gene_id_version, ensembl_gene_id, ensembl_transcript_id_version and external_gene_name and these *MUST* be named as is
# If the gtf has the correct format then provide set args0 to 'nothing' and provide path to gtf in args1 to subsset the required columns. Set args2 to 'NULL'
# If no gtf, set args1 to 'download', args1 to 'NULL' and provide biomart dataset and ENSEMBL version in args2, such that,
# .. element 1 is named 'ds' (biomart dataset) and element 2 is named 'version' (ensembl version) e.g c('ds'='hsapiens_gene_ensembl','version'='92'). Argument restricted to 2 elements
# If gtf needs formatting set args0 to 'format' and provide path to get in args1. Set args2 such that,
# .. if any colum is missing then element name is 'missing' and element value is the name of the missing column (as in the definition above)
# .. .. In addition element 'ds' and element 'version' need to point to biomart dataset and ensembl version to be use for retrieving the missing information.
# .. .. Add one 'missing' element per missing column.
# .. .. add the column anchor (element named 'anchor') for merging both the provided gtf and the retrieved information
# .. .. .. it performs a left join on the provided gtf, that is, any ID that *is not* in the gtf is NOT kept
# .. If any columns in the provided gtf needs renaming to fullfil the above definition, add one element per field that needs renaming, such that
# .. .. element name is the original name present in the provided gtf and element value is the new name as defined above
# .. Argument can have more than 2 elements


  if(!"package:tidyverse" %in% search()) library("tidyverse")
  if(!"package:rtracklayer"%in%search()) library("rtracklayer")
  if(!"AnnotationHub"%in%search()) library("AnnotationHub")
  if(!"package:biomaRt"%in%search()) library("biomaRt")

  if(args5=='download'){
    cat("..downloading data\n")
    gid=getAnn_f(dataset=args2['ds'], version=args2['version'])

    if(args3=="yes") {
      saveRDS(gid,file=paste(args4,".rds",sep=""))
      cat("..saving",paste(args4,".rds",sep=""),"\n")
    }
  }

  if(args5=='nothing'){
    gtf= readGFF(args1,columns = GFFcolnames())
    gid=gtf %>% select ("ensembl_gene_id_version","ensembl_gene_id","ensembl_transcript_id_version","external_gene_name")
    gid=unique(gid)
  }

  if(args5=='format'){
    cat("..formatting GTF\n")
    gtf= readGFF(args1,columns = GFFcolnames())

    for(i in which(!names(args2)%in%c("missing","ds","version","anchor"))){
        colnames(gtf)[colnames(gtf)==names(args2)[i]]=args2[i]
    }

    if(any(names(args2)=="missing")) {
      cat("..downloading data\n")
      gid=getAnn_f(dataset=args2['ds'], version=args2['version'])
    }

    idx=colnames(gtf)%in%c("ensembl_gene_id_version","ensembl_gene_id","external_gene_name","ensembl_transcript_id_version")
    tmp=gtf[,idx] %>% as_tibble(); colnames(tmp)=colnames(gtf)[idx]; rm(idx)
    tmp=unique(tmp)

    if(exists("gid")) {gid=left_join(tmp, gid, by=paste(args2["anchor"])); comment(tmp)="suffix 'x' refers to provided gtf, suffix 'y' refers to downloaded from biomart data"}

    if(args3=="yes"){
      saveRDS(gid, file=paste(args4,".rds",sep=""))
      cat(".. saving",paste(args4,".rds",sep=""),"\n")
    }
  }
  return(gid)
}

rename_rows_seuratF <- function(obj,index){
	# obj is a seurat object with AT LEAST a count matrix in the RNA assay that NEEDS to be named 'counts'
	# index is a dataframe [external_gene_name, ensembl_gene_id]
	# the mapping is directional external_gene_name --> ensembls_gene_id
	if(!"package:dplyr" %in% search()) library("dplyr")
	if(!"package:tidyr" %in% search()) library("tidyr")
	
	# remove genes with no ensembl information
	index=index[match(rownames(obj),index$external_gene_name),c("external_gene_name","ensembl_gene_id")] %>% as_tibble() %>% drop_na()
	idx1=match(index$external_gene_name,rownames(obj))
	cat(".. check 1\n")
	
	obj=obj[idx1,]
	
	ctm=GetAssayData(object=obj,assay="RNA",slot="counts")
	rownames(ctm)=index$ensembl_gene_id
	rm(idx1)
	cat(".. counts\n")
	
	tmp = CreateSeuratObject(counts=ctm,assay="RNA")
	cat("..seurat object\n")
	
	assay=names(obj)
	
	for(i in assay){
		asstp=class(obj[[paste(i)]])
		if(asstp=="Assay"){
			cat(".. assays\n")
			
			# normalized
			datanm=GetAssayData(object=obj,assay=i,slot="data")
			rownames(datanm)=index$ensembl_gene_id[match(rownames(datanm),index$external_gene_name)]
			# scaled
			snm=GetAssayData(object=obj,assay=i,slot="scale.data")
			rownames(snm)=index$ensembl_gene_id[match(rownames(snm),index$external_gene_name)]
			
			if(i=="RNA"){
				tmp = SetAssayData(tmp, assay=i, slot="data", new.data=datanm)
				tmp = SetAssayData(tmp, assay=i, slot="scale.data", new.data=snm)
			} else {
				nA=CreateAssayObject(data=datanm)
				tmp[[paste(i)]]=nA
				tmp = SetAssayData(tmp, assay=i, slot="scale.data",new.data=snm)
			}
			rm(datanm); rm(snm)
			
			vF=VariableFeatures(obj,assay=i)
			vF=index$ensembl_gene_id[match(vF,index$external_gene_name)]
			VariableFeatures(tmp,assay=i)=vF
			rm(vF)
		}
		if(asstp=="Graph"){
			cat(".. graphs\n")
			tmp[[paste(i)]]=obj[[paste(i)]]
		}
		if(asstp=="DimReduc"){
			cat(".. reduced dims\n")
			fL=Loadings(obj,reduction=i)
			rownames(fL)=index$ensembl_gene_id[match(rownames(fL),index$external_gene_name)]
			tmp[[paste(i)]]=CreateDimReducObject(embeddings=Embeddings(obj,reduction=i),
				loadings=fL,
				stdev = obj[[paste(i)]]@stdev,
				jackstraw = obj[[paste(i)]]@jackstraw,
				misc = obj[[paste(i)]]@misc,
				key=paste(i,"_",sep=""),
				assay=obj[[paste(i)]]@assay.used)
		}
	}
	if(ncol(obj@meta.data)>0) tmp = AddMetaData(tmp, obj@meta.data)
	#tmp@neighbors=obj@neighbors
	tmp@version=obj@version
	tmp@commands=obj@commands
	tmp@tools=obj@tools
	DefaultAssay(tmp)=DefaultAssay(obj)
	Idents(tmp)=Idents(obj)
	return(tmp)
}
#source("/castor/project/home/mararc/bin/as.Seurat.SingleCellExperiment_develBranch_commit723b066_Dec7-2022.R")

rename_rows_seurat_v5F <- function(obj,index){
	# obj is a seurat object with AT LEAST a count or data matrix in the RNA assay that NEEDS to be named 'counts'
	# index is a dataframe [external_gene_name, ensembl_gene_id]
	# the mapping is directional external_gene_name --> ensembls_gene_id
	if(!"package:dplyr" %in% search()) library("dplyr")
	if(!"package:tidyr" %in% search()) library("tidyr")
	
	# remove genes with no ensembl information
	index=index[match(rownames(obj),index$external_gene_name),c("external_gene_name","ensembl_gene_id")] %>% as_tibble() %>% drop_na()
	idx1=match(index$external_gene_name,rownames(obj))
	cat(".. check 1\n")
	
	obj=obj[idx1,]
	rm(idx1)
	
	# counts
	if(obj[["RNA"]]$counts %>% dim() %>% sum()>0) {
		ctm=GetAssayData(object=obj,assay="RNA",layer="counts")
		rownames(ctm)=index$ensembl_gene_id
		cat(".. RNA counts\n")
	} else {ctm=NULL}
	
	# normalized
	if(obj[["RNA"]]$data %>% dim() %>% sum()>0) {
		datanm=GetAssayData(object=obj,assay="RNA",layer="data")
		rownames(datanm)=index$ensembl_gene_id[match(rownames(datanm),index$external_gene_name)]
		cat(".. RNA data\n")
	} else {datanm = NULL}
	
	# scaled
	if(obj[["RNA"]]$scale.data %>% dim() %>% sum()>0) {
		snm=GetAssayData(object=obj,assay="RNA",layer="scale.data")
		rownames(snm)=index$ensembl_gene_id[match(rownames(snm),index$external_gene_name)]
		cat(".. RNA scale.data\n")
	} else {snm = NULL}
	
	# build object
	if(sum(dim(ctm))>0) { 
		tmp = CreateSeuratObject(counts=ctm,assay="RNA")
	} else {
		tmp = CreateSeuratObject(counts=datanm,assay="RNA")
		cat(".. no count matrix, using data as a placeholder\n")
		}
	if(sum(dim(datanm))>0) tmp = SetAssayData(tmp, assay="RNA", layer="data", new.data=datanm)
	if(sum(dim(snm))>0) tmp = SetAssayData(tmp, assay="RNA", layer="scale.data", new.data=snm)
	rm(ctm); rm(datanm); rm(snm); gc()
	cat("..seurat object\n")
	
	vF=VariableFeatures(obj,assay="RNA")
	vF=index$ensembl_gene_id[match(vF,index$external_gene_name)]
	VariableFeatures(tmp,assay="RNA")=vF
	rm(vF)
	
	assay=names(obj)[!names(obj)%in%"RNA"]
	for(i in assay){
		asstp=class(obj[[paste(i)]])
		if(asstp=="Assay"){
			cat(".. assays\n")
			# counts
			cts=GetAssayData(object=obj,assay=i,layer="counts")
			rownames(cts)=index$ensembl_gene_id[match(rownames(cts),index$external_gene_name)]
			# normalized
			datanm=GetAssayData(object=obj,assay=i,layer="data")
			rownames(datanm)=index$ensembl_gene_id[match(rownames(datanm),index$external_gene_name)]
			# scaled
			snm=GetAssayData(object=obj,assay=i,layer="scale.data")
			rownames(snm)=index$ensembl_gene_id[match(rownames(snm),index$external_gene_name)]
			
			nA=CreateAssayObject(counts=cts)
			tmp[[paste(i)]]=nA 
			if(sum(dim(datanm))>0) tmp = SetAssayData(tmp, assay=i, layer="data",new.data= datanm)
			if(sum(dim(snm))>0) tmp = SetAssayData(tmp, assay=i, layer="scale.data",new.data=snm)
			
			rm(datanm); rm(snm)
			
			vF=VariableFeatures(obj,assay=i)
			if(length(vF)>0) {
				vF=index$ensembl_gene_id[match(vF,index$external_gene_name)]
				VariableFeatures(tmp,assay=i)=vF
				rm(vF)
			}
		}
		if(asstp=="Graph"){
			cat(".. graphs\n")
			tmp[[paste(i)]]=obj[[paste(i)]]
		}
		if(asstp=="DimReduc"){
			cat(".. reduced dims\n")
			fL=Loadings(obj,reduction=i)
			rownames(fL)=index$ensembl_gene_id[match(rownames(fL),index$external_gene_name)]
			tmp[[paste(i)]]=CreateDimReducObject(embeddings=Embeddings(obj,reduction=i),
				loadings=fL,
				stdev = obj[[paste(i)]]@stdev,
				jackstraw = obj[[paste(i)]]@jackstraw,
				misc = obj[[paste(i)]]@misc,
				key=paste(i,"_",sep=""),
				assay=obj[[paste(i)]]@assay.used)
		}
	}
	
	if(length(obj@misc)>0) tmp@misc = obj@misc
	if(ncol(obj@meta.data)>0) tmp = AddMetaData(tmp, obj@meta.data)
	
	#tmp@neighbors=obj@neighbors
	tmp@version=obj@version
	tmp@commands=obj@commands
	tmp@tools=obj@tools
	DefaultAssay(tmp)=DefaultAssay(obj)
	Idents(tmp)=Idents(obj)
	return(tmp)
}

my_readGTF <- function(gtf){
	# input (gtf) is the path to a GTF annotation file
	# rtracklayer::import.gff() reads a gtf. But, the image doesn't have that package
	# outputs and index which gene names have been adapted to seurat format
	# Field 9 in gtf has to have "gene_id" and "unique_gene_name"
	# # If "unique_gene_name" is not present, adapt the code to read "gene_name" instead 
	gtf=read.table(file=gtf,header= FALSE, sep="\t")
	gtf=gtf[gtf[,3]=="exon",]
	# Get ensembl_gene_id and unique_gene_name - don't need to check for ENSEMBL ambiguity because I am mapping based on unique_gene_name
	gtfidx=lapply(gtf[,9], function(x) {strsplit(x,split=";")[[1]]%>%.[grep(x=.,pattern="unique_gene_name|gene_id")]}) %>% do.call(rbind,.) %>% trimws(.) %>% data.frame(.)
	colnames(gtfidx)=c("ensembl_gene_id_version","external_gene_name")
	gtfidx$ensembl_gene_id_version=sapply(gtfidx$ensembl_gene_id_version, function(x) {strsplit(x,split="gene_id ")[[1]][2] %>% trimws(.)})
	gtfidx$external_gene_name=sapply(gtfidx$external_gene_name, function(x) {strsplit(x,split="unique_gene_name ")[[1]][2]%>% trimws(.)})
	gtfidx=data.frame(unique(gtfidx))
	gtfidx$ensembl_gene_id=sapply(gtfidx$ensembl_gene_id_version, function(x) strsplit(x,split=".",fixed=TRUE)[[1]][1])
	# Because Seurat replaces "_" with "-", names are different from GTF. I need to modify names to fit as needed
	gtfidx$external_gene_name=gsub("_","-",gtfidx$external_gene_name)
	# Remove ambiguities
	tmp=gtfidx%>% group_by(ensembl_gene_id)%>%summarize(n=n_distinct(external_gene_name)) %>% dplyr::filter(n>1) %>% dplyr::select(ensembl_gene_id) %>% unlist()
	if(length(tmp)>0) gid = gid %>% dplyr::filter(! ensembl_gene_id %in% tmp)
	rm(tmp)
	tmp=gtfidx%>% group_by(external_gene_name)%>%summarize(n=n_distinct(ensembl_gene_id)) %>% dplyr::filter(n>1) %>% dplyr::select(external_gene_name) %>% unlist()
	if(length(tmp)>0) gid = gid %>% dplyr::filter(! external_gene_name %in% tmp)
	rm(tmp)
	rm(gtf)
	return(gtfidx)
}

readFileF <- function(input_file){
	# Input_file can be binary (RDS, RData) or a flat file ("txt","tsv","csv","bed")
	# If RData the file should have the objects names as "sce" for SingleCellExperiment or "seudat" for Seurat
	# If flat file, they must have a header
	# It can also have extra gene metadata stored in "rowdata" and miscelaneous in "misc"
	
	if(!"package:dplyr" %in% search()) library("dplyr")
	if(!"package:tidyr" %in% search()) library("tidyr")

	fltype=strsplit(input_file,split=".",fixed=TRUE)[[1]]%>% .[length(.)] %>% tolower(.) %>% unlist()
	
	if(fltype=="rds"){ fileBCK=readRDS(input_file) }
	if(fltype=="rdata") {
		attach(input_file)
		if(exists("sce")) {fileBCK=sce}
		if(exists("seudat")) {fileBCK=seudat}
		if(exists("rowdata")) {rowdataBCK=rowdata}
		if(exists("misc")) {miscBCK=misc}
		detach()
	}
	if(fltype%in%c("txt","tsv","csv","bed")){ fileBCK=read.table(input_file,header=TRUE) }
	if(fltype=="gtf"){fileBCK=my_readGTF(input_file)}
	if(!fltype%in%c("txt","tsv","csv","bed","gtf","rds","rdata")){ fileBCK=eval(parse(text=input_file)) }
	
	if(exists("fileBCK")&exists("rowdataBCK")&exists("miscBCK")) tmp=list(dat=fileBCK,rowdata=rowdataBCK,misc=miscBCK)
	if(exists("fileBCK")&exists("rowdataBCK")&!exists("miscBCK")) tmp=list(dat=fileBCK,rowdata=rowdataBCK)
	if(exists("fileBCK")&exists("miscBCK")&!exists("rowdataBCK")) tmp=list(dat=fileBCK,misc=miscBCK)
	if(exists("fileBCK")&!exists("rowdataBCK")&!exists("miscBCK")) tmp=list(dat=fileBCK)
	
	return(tmp)
}

contingency_tablesF <- function(data,grouping_variable,testing_variable){
	# data is a data.frame with the two variables needed to create the contingency tables 
	# grouping_variable is the column name in data where the variable of interest is
	# testing variable is the column name in data where the variable to be tested is
	
	if(!"package:dplyr" %in% search()) library("dplyr")
	
	DD=data
	GV=grouping_variable
	TV=testing_variable
	
	idx=expand.grid(DD[,GV]%>%unlist()%>%unique(),DD[,TV]%>%unlist()%>%unique()) %>%unique()
	idx$Var1=as.character(idx$Var1); idx$Var2=as.character(idx$Var2)
	tmp=lapply(1:nrow(idx), function(i) {
		tmp=table(DD[,GV]%>%unlist(),DD[,TV]%>%unlist()) %>% as.data.frame.matrix() %>% mutate(REST=rowSums(.[,colnames(.)!=unlist(idx[i,2]),drop=FALSE]))%>% select(unlist(idx[i,2]),REST) 
		tmp1=tmp %>% filter(!rownames(.)%in%unlist(idx[i,1]))%>%colSums()
		tmp= bind_rows(tmp %>% filter(rownames(.)==unlist(idx[i,1])),tmp1) 
		rm(tmp1)
		colnames(tmp)=c(unlist(idx[i,2]),paste("not-",unlist(idx[i,2]),sep=""))
		rownames(tmp)=c(unlist(idx[i,1]),"other")
		return(tmp)
	})
	names(tmp)=paste(GV,"-",unlist(idx[,1]),"_",TV,"-",unlist(idx[,2]),sep="")
	return(tmp)
}

choose_qual_color_f <- function(n, seed=1234){
	# n is an integer of the number of colors to be returned
	# needs RColorBrewer
	
	if(!"package:RColorBrewer" %in% search()) library("RColorBrewer")
	
	col_pals=brewer.pal.info[brewer.pal.info$category=='qual',]
	col_vec=unlist(mapply(brewer.pal, col_pals$maxcolors, rownames(col_pals)))
	
	set.seed(seed)
	col_vec=sample(col_vec,n)
	return(col_vec)
}

add_tpmF <- function(obj,glidx){
	# obj is either a Seurat or SCE object
	# glidx a named vector of length equal to nrow(obj) sorted as obj, such that all(names(glidx)==rownames(obj))
	# glidx contains gene lengths
	# NOTE: because these are sparse matries, the print return by seurat would be altered, but the matrices sizes remained
	
	if(!"package:Seurat" %in% search()) library("Seurat")
	if(!"package:SingleCellExperiment" %in% search()) library("SingleCellExperiment")
	
	if(class(obj)=="Seurat") {
		mocksce=as.SingleCellExperiment(obj)
		tpm=scater::calculateTPM(mocksce,assay.type="counts",lengths=glidx)
		obj[["TPM"]]=CreateAssayObject(data=tpm)
		rm(mocksce); rm(tpm)
	} 
	if(class(obj)=="SingleCellExperiment") {
		assay(obj,"TPM")=scater::calculateTPM(obj,assay.type="counts",lengths=glidx)
	}
	return(obj)
}

longTranscript <- function(dat) {
	## NEED dplyr and scde
	## Method to retain only the longest spanning trancript
	## INPUT is a data frame with columns named: gene, chr, start, end
	gene_loc<- dat
	gene_loc<- distinct(gene_loc)
	## For duplicated transcripts -by chromosome-, keep only the transcript with the smaller start position and the transcript with the largest end position. 
	## When several transcripts have the same start or same end, keep the longest one
	## These are the  transcripts that span the extremes, that is, the longest interval 
	dup<-  as.data.frame( group_by(gene_loc,gene,chr) %>% filter(start == min(start) | end == max(end)) )
	dup <- dup %>% mutate(dist = end - start)
	## When more than 1 transcript have the same smaller starting position, keep longest transcript
	maxDS <- data.frame( summarise(group_by(dup,chr,gene,start),dist = max(dist)) ) 
	dup <- semi_join(dup, maxDS, by=c("chr", "gene", "start", "dist"))
	## When more than 1 transcript have the same largest ending position, keep the longest transcript
	maxDE<- data.frame( summarise(group_by(dup,chr,gene,end), dist = max(dist)) ) 
	dup <- semi_join(dup, maxDE, by=c("chr", "gene", "end", "dist")) 
	dup<- dup[,c("gene","chr","start","end")]
	## split table by chromosome
	chr_l<- lapply(unique(dup$chr), function(x) dup[dup$chr == x,]) 

	newT<- vector("list", length(chr_l))
	for(i in 1:length(chr_l)) {
		## Get transcripts that are duplicated
		## TRANSCRIPTS OVERLAP: If start position of duplicated transcript is smaller than end position of trancript above, set new start position to start position of gene above
		## TRANSCRIPTS DON'T OVERLAP: if start position of duplicated trancripts is larger than end position of trancript above, set new start position to old start position
		dat<- chr_l[[i]]
		dat <- arrange(dat, gene, start) # sort based on gene and starting position
		idx1<- which(duplicated(dat$gene))
		idx2<- which(!duplicated(dat$gene))
		dat[idx1,"newS"]<-ifelse(dat[idx1,"start"] < dat[idx1-1,"end"], dat[idx1-1, "start"], dat[idx1,"start"])
		dat[idx2,"newS"]<- dat[idx2,"start"]
		dat[idx1,"newE"]<-dat[idx1,"end"]
		dat[idx2,"newE"]<- dat[idx2,"end"]
		dat<- dat[,c("gene","chr","newS","newE")]
		colnames(dat)<- c("gene","chr","start","end")
		newT[[i]]<- dat 
							}
							
	newT<- do.call(rbind,newT)	
	newT<- newT %>% mutate(dist= end-start)	
	maxD <- data.frame(summarise(group_by(newT,gene), dist = max(dist)))	## Get only the longest trancript
	newT<- semi_join(newT, maxD, by=c("gene","dist"))
	## Ignore trancripts that have the exact same longitude but are in different locations
	idx3<-newT$gene[duplicated(newT$gene)]
	newT<- newT[!(newT$gene %in% idx3),-5]
	return(newT)
}

getCNVmat<- function(WD, cellgroupings, stateprob, cnvreg, genesincnv ) {
	
	# simplified version of getCNVdat5(). Created in 2023-07-17. Does not take information on chromosome arms (p and q) or clusters
	# Function to retrieve the data files where the CNV information is, that are generated by an infercnv::run()
	# WD is the directory passed to infercnv::run(out_dir= ) where the files will be written
	# Needs dplyr
	# Needs a number of files to be input as strings:
	# .. stateprob: infercnv CNV_State_Probabilities.dat file
	# .. cnvreg: infercnv .rpred_cnv_regions.dat
	# .. genesincnv: infercnv .pred_cnv_genes.dat
	# .. cellgroupings: infercnv .cell_groupings file
	#---------------------------------------------------------------------------------------------#
	
	if(!"package:dplyr" %in% search()) library("dplyr")
	#---------------------------------------------------------------------------------------------#
	
	stateprob<-stateprob # paste(WD,"/BayesNetOutput.HMMi6.leiden.hmm_mode-subclusters/CNV_State_Probabilities.dat", sep="")
	cnvreg<- cnvreg # paste(WD,"HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.1.pred_cnv_regions.dat",sep="") 
	genesincnv<- genesincnv #paste(WD,"HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.1.pred_cnv_genes.dat",sep="")
	cellgroupings<- cellgroupings #paste(WD,"/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings",sep="")
	#---------------------------------------------------------------------------------------------#
	
	## Get the probability of each state for each CNV region
	cat("Getting state probabilities for significant CNVs:\n",stateprob,"\n")
	regprob<- read.table( stateprob, header = TRUE)
	fregprob<- data.frame(Prob= apply(regprob,2,function(x) max(x)), State= apply(regprob,2, function(x) names(x)[x == max(x)])) ## format regprob
	fregprob<- data.frame(fregprob, cnv_name= gsub(".","-",rownames(fregprob), fixed=TRUE))
	fregprob$State<-as.numeric(sapply(lapply(fregprob$State, function(x) strsplit(as.character(x),":")),function(z)  z[[1]][2]))

	## Get the cnv regions
	cat("Getting the CNV regions:\n",cnvreg,"\n")
	regpred<- read.table(cnvreg, header=TRUE)
	regpred$cell_group_name<- as.character(regpred$cell_group_name)
	
	## Get the cell groupings
	cat("Getting the cell groupings:\n",cellgroupings,"\n")
	cellgroup<- read.table(cellgroupings, header=TRUE)
	cellgroup$cell_group_name<- as.character(cellgroup$cell_group_name)

	## Get genes in CNV
	cat("Getting genes within each CNV:\n",genesincnv,"\n")	
	genes<-read.table(genesincnv, header=TRUE)

	# cnv coordinates
	mer<-left_join(regpred, fregprob, by="cnv_name")
	mer<-right_join(mer, cellgroup, by="cell_group_name", relationship="many-to-many") # a right join keeps all cells. Wherever there are NAs that means the cell did not have significant CNVs
	mer$Prob[is.na(mer$cnv_name)]<-0
	mer$State[is.na(mer$cnv_name)]<-3
	mer$cnv_name[is.na(mer$cnv_name)]<-"none"

	#mer<-left_join(mer,clus,by="cell")
	mer<-mer[,!colnames(mer)%in%c("state")]; colnames(mer)[which(colnames(mer)=="State")]="state"
	
	# gene coordinates
	mer1<-left_join(genes, fregprob, by=c("gene_region_name"="cnv_name"))
	mer1<-right_join(mer1, cellgroup, by="cell_group_name", relationship="many-to-many") # a right join keeps all cells. Wherever there are NAs that means the cell did not have significant CNVs
	#mer1$Prob[is.na(mer1$gene_region_name)]<-0
	#mer1$State[is.na(mer1$gene_region_name)]<-3
	#mer1$gene_region_name[is.na(mer1$gene_region_name)]<-"none"
	#mer1$gene[is.na(mer1$gene_region_name)]<-"none"
	mer1<-mer1[,!colnames(mer1)%in%c("state")]; colnames(mer1)[which(colnames(mer1)=="State")]="state"

	return(list(cnvcoord=mer, genecoord=mer1))
}

gene_set_score <- function(genes, obj){
	# provide a vector with gene ids (need to be in rownames(obj))
	# seurat object
	
	Gidx=genes
	Sobj=obj
	
	# Get per cell mean expression
	mean_exp=colMeans(x = Sobj[['RNA']]@data[Gidx,],na.rm=TRUE)
	
	# Store scores
	if(all(names(x = mean_exp) == rownames(Sobj@meta.data))) {
			cat("adding per cell signature score to metadata\n")
			Sobj@meta.data$gene_set_score=mean_exp
	}
	return(Sobj)
}

getRdata <- function(rdatafile, robject) {
	E = new.env()
	load(rdatafile, envir=E)
	
	if(length(robject)>1){ 
		return(lapply(robject, function(x) get(x, envir=E, inherits=FALSE)))
	} else {
		return(get(robject,envir=E, inherits=FALSE))
	}
} 

maxs_f=function(row,n){
	# identify n maximum values of a vector an get name
	sortv=sort(unlist(row),decreasing=TRUE)
	maxs=sortv[1:n]
	nms=sort(names(maxs))
	nms=paste(nms,collapse="_")
	return(nms)
}

max_part_distance=function(row,max_group_size){
	# Calculate the distance between groups
	calculate_distance=function(group1, group2){
	mean(group1) - mean(group2)
	#abs(sum(group1) - sum(group2))
	#sqrt((mean(group1)-mean(group2))^2) euclidean make sense when n-dimensional (vectors) but it is the same as absolute distance in the 2 dimensional space
	}

	# get all possible partitions
	n=max_group_size  #length(row)
	parts=list()
	for(i in 1:n){
		combs=combn(row,i,simplify=FALSE)
		parts=c(parts,combs)
	}

	# calculate distance for each partition
	dist_list = lapply(parts, function(x) {
		g1 = x %>% unlist()
		g2 = row[!names(row)%in%names(g1)] %>% unlist()   #setdiff(row %>% unlist(),x)
		R=list(group1=g1, group2=g2, abs_dist=abs(calculate_distance(g1,g2)))
		return(R)
		}
	)

	# find partition with max dist
	distances=sapply(dist_list,function(x) x$abs_dist)
	max_distance=which.max(distances)
	
	# get group with max mean score
	grp=dist_list[[max_distance]]
	m1=mean( as.vector(as.numeric(grp[["group1"]])) )
	m2=mean( as.vector(as.numeric(grp[["group2"]])) )
	if(which.max(c(m1,m2))==1) {G=list(group=grp[["group1"]],label=paste(names(grp[["group1"]])%>%sort(),collapse="_"),mean_score=m1)}
	if(which.max(c(m1,m2))==2) {G=list(group=grp[["group2"]],label=paste(names(grp[["group2"]])%>%sort(),collapse="_"),mean_score=m2)}
	
	return(G)
}
max_part_distance1=function(row,max_group_size){ # only diff from above is collapsing with "**"
	# Calculate the distance between groups
	calculate_distance=function(group1, group2){
	mean(group1) - mean(group2)
	#abs(sum(group1) - sum(group2))
	#sqrt((mean(group1)-mean(group2))^2) euclidean make sense when n-dimensional (vectors) but it is the same as absolute distance in the 2 dimensional space
	}

	# get all possible partitions
	n=max_group_size  #length(row)
	parts=list()
	for(i in 1:n){
		if(length(row)>i){
			combs=combn(row,i,simplify=FALSE)
			parts=c(parts,combs)
		}
	}

	# calculate distance for each partition
	dist_list = lapply(parts, function(x) {
		g1 = x %>% unlist()
		g2 = row[!names(row)%in%names(g1)] %>% unlist()   #setdiff(row %>% unlist(),x)
		R=list(group1=g1, group2=g2, abs_dist=abs(calculate_distance(g1,g2)))
		return(R)
		}
	)

	# find partition with max dist
	distances=sapply(dist_list,function(x) x$abs_dist)
	max_distance=which.max(distances)
	
	# get group with max mean score
	grp=dist_list[[max_distance]]
	m1=mean( as.vector(as.numeric(grp[["group1"]])) )
	m2=mean( as.vector(as.numeric(grp[["group2"]])) )
	if(which.max(c(m1,m2))==1) {G=list(group=grp[["group1"]],label=paste(names(grp[["group1"]])%>%sort(),collapse="**"),mean_score=m1)}
	if(which.max(c(m1,m2))==2) {G=list(group=grp[["group2"]],label=paste(names(grp[["group2"]])%>%sort(),collapse="**"),mean_score=m2)}
	
	return(G)
}

predictF=function(X,y,predictor_sets,alpha=1,n_splits){
	
}

sample_geneset=function(vec, num_samples, sample_size){
	# function to sample genes at least once while keeping duplicates to a minimum
	chunk_size=length(vec)%/%num_samples
	remainder=length(vec)%%num_samples
	sampled_sets=list()
	sampled_items=integer(0)
	
	for(i in 1:num_samples){
		start_idx=(i-1)*chunk_size + 1
		end_idx=start_idx + chunk_size -1
		
		if(i<=remainder){
			end_idx=end_idx+1
		}
		
		chunk=vec[start_idx:end_idx]
		samp=sample(chunk,sample_size)
		sampled_sets[[i]]=abs(samp)
		#sampled_items=c(sampled_items,samp)
	}
	return(sampled_sets=sampled_sets)
}

fit_modelF=function(
	X,
	y,
	test_size=0.2,
	nfolds=10, # for lambda default grid
	type.measure="deviance",
	q_stabsel=20
	){
	# this function takes a set of predictors and trains a lasso model
	
	# Split into training and testing
	train_idx=createDataPartition(y,p=1-test_size,list=FALSE)
	X_train=X[train_idx,]
	y_train=y[train_idx]
	X_test=X[-train_idx,]
	y_test=y[-train_idx]

	# Fit model
	cv=cv.glmnet(
		x=X_train,
		y=y_train,
		alpha=1, # 0= ridge 0.5=elastic net 1=lasso
		family="binomial",
		nfolds=nfolds,
		type.measure=type.measure,
		parallel=TRUE
	)
	best_lambda=cv$lambda.min
	y_probs=predict(cv,newx=X_test,type="response",s=best_lambda)
	y_pred=ifelse(y_probs>0.5,1,0)
	y_pred=factor(y_pred,levels=levels(y_test))
	raw_coeffs=coef(cv,s=best_lambda)

	# Stability selection
	if(FALSE){
	m1=function(x,y,q,...){
		cv_fit=cv.glmnet(x,y,alpha=1,family="binomial")
		return(list(as.vector(coef(cv_fit,s="lambda.min"))))
	}
	stasel=stabsel(x=X_train,y=y_train,fitfun=m1,cutoff=0.75,q=q_stabsel)
	}
	
	
	# Confusion matrix
	cm=confusionMatrix(y_pred,y_test)
	
	# More metrics
	sens=cm$byClass["Sensitivity"]
	spec=cm$byClass["Specificity"]
	ppv=cm$byClass["Pos Pred Value"]
	nppv=cm$byClass["Neg Pred Value"]
	accur=cm$overall["Accuracy"]
	f1=cm$byClass["F1"]
	auc=roc(y_test,as.numeric(y_pred))$auc # control = reference level = baseline against which the positive is compared

	return(
		list(
			model=cv,
			raw_preds=y_pred,
			sensitivity=sens,
			specificity=spec,
			ppv=ppv,
			accuracy=accur,
			f1=f1,
			auc=auc,
			coefficients=coeffs
		)
	)
}

glmnet.lasso_logistic=function (x, y, q, type = c("conservative", "anticonservative"), ...) {
	# function to fit a lasso logistic regression.
	# inspired by stabsel::lars.lasso
    if (!requireNamespace("glmnet", quietly = TRUE)) 
        stop("Package ", sQuote("glmnet"), " needed but not available")
    if (is.data.frame(x)) {
        message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
        x = model.matrix(~. - 1, x)
    }
    if ("lambda" %in% names(list(...))) 
        stop("It is not permitted to specify the penalty parameter ", 
            sQuote("lambda"), " for lasso when used with stability selection.")
    type = match.arg(type)
    if (type == "conservative") 
        fit = suppressWarnings(glmnet::glmnet(x, y, pmax = q,
			family="binomial", alpha=1, nfolds=10, parallel=TRUE...))
    if (type == "anticonservative") 
        fit = glmnet::glmnet(x, y, dfmax = q - 1,
			family="binomial", alpha=1, nfolds=10, parallel=TRUE...)
    selected = predict(fit, type = "nonzero")
    selected = selected[[length(selected)]]
    ret = logical(ncol(x))
    ret[selected] = TRUE
    names(ret) = colnames(x)
    cf = fit$beta
    sequence = as.matrix(cf != 0)
    return(list(selected = ret, path = sequence))
}

fastDotPlot=function(obj,parm){
	args41=lapply(1:length(parm), function(i) {
		x=parm[i] %>% unlist()
		tmp=x[x%in%rownames(obj)]
		names(tmp)=NULL
		return(tmp)
	}
	); names(args41) = names(parm)
	args41= args41[!sapply(args41,function(x) length(x)==0)]

	tmp=DotPlot(obj,assay="RNA",
		features=args41[1:length(args41)],
		group.by="custom_cluster_counts", 
		cols=c("blue","red"),
		dot.scale=10,
		scale=TRUE)+ #
		theme(
			axis.text.x=element_text(
				angle=45,
				hjust=1,
			#	size=6
				),
			axis.text.y=element_text(
			#	size=6,
				hjust=1,
				angle=25
				),
			#strip.text=element_text(size=6),
			strip.background=element_rect(
				fill="lightgrey",
				color="black"
				),
			#legend.text=element_text(size=6),
			#legend.title=element_text(size=8),
			panel.spacing=unit(0.1,"lines"),
			axis.title.x=element_blank(),
			panel.border=element_rect(color="grey",fill=NA,linewidth=1),
			legend.position="bottom",
			legend.direction="horizontal",
			legend.box="horizontal",
			plot.margin=unit(c(1,1,1,1),"cm")
		)+
		guides(
			color=guide_colourbar(title="mean scaled\nexpression"), # <- or normalized when scale=FALSE
			size=guide_legend(title="percent\nexpressed",nrow=2)
		)+
		ylab("tumor clusters")
		
		return(tmp)
}
