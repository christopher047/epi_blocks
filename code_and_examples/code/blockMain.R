
#### parse bismark data  ####
source("parseFromBismark.R") 
n          <- 4 	# n minimum reads, less than n than read is set NA for impution
n_imp      <- 70 	#n_imp minimum percentage of sites with coverage values greater than n
files      <- list.files(path="../bismark_data/", pattern="*.gz$", full.name=T)
parseBismark(files, n=4, n_imp=70, output_dir="../chrs/")

### results ###
# look at created files files ../chrs/chr2_tot_meth.tab and ../chrs/chr2_tot_cov.tab 
# rownames are chromosome then position seperated by a period or chr2.10000 
# the rest are values 

### paramaters ####
achr="chr2"  		#chromosome to work on 
input_dir="../chrs/"    #directory where chr2_tot_meth.tab and chr2_tot_cov.tab files are located
min.seg=20		#minimum number of CpG sites in a segment 
min.block=10	        #minimum number of CpG sites in a block 
bwd=300			#bandwidth decay paramater see blockMethods methylDist
hclust=0.4		#hierarchal clustering cut tree cutoff see blockMethods endreCluster
iqr_cutoff=10	        #filter Range (max -min) divided by IQR (75 percentile - 25 percentile) less than 10
nb=5			#number of neighbours for KNN to use in imputation 	
ncores=16		#number of cores to use 

### run for all samples  ####
source("blockClassDefinition.R") 
my_start       <- Sys.time() 
chr2           <- runChromosome(achr=achr, input_dir=input_dir,  min.seg=min.seg, bwd=bwd, min.block=min.block, 
				hclust=hclust, iqr_cutoff=iqr_cutoff, nb=nb, ncores=ncores)

cat(paste("finished", Sys.time()-my_start,"\n"), sep="\n") 


### results ####
imp_meth <- chr2$imp_meth	#imputed methylation matrix 
blocks   <- chr2$blocks		#blocks Granges object
raw_cov  <- chr2$raw_cov	#raw coverage matrix to check 
index    <- chr2$index		#index to link block name to imputed methylation matrix rownames 
segs     <- chr2$segs		#list of segments 

### save ###
write.table(as.data.frame(imp_meth), "../results/chr2_imp_meth.tab", sep="\t") 
write.table(as.data.frame(blocks), "../results/chr2_blocks.tab", sep="\t") 
write.table(as.data.frame(index), "../results/chr2_index.tab", sep="\t") 

### load ###
index    <- GRanges(read.table("../results/chr2_index.tab", row.names=1, header=T, sep="\t"))
blocks   <- GRanges(read.table("../results/chr2_blocks.tab", row.names=1, header=T, sep="\t"))
imp_meth <- read.table("../results/chr2_imp_meth.tab", row.names=1, header=T, sep="\t")

### overlap #####
source("getCpgIslands.R") 
cpg      <- getCpgIslands() 
ov       <- findOverlaps(blocks, cpg)
pct_cpg  <- length(unique(queryHits(ov)))/length(blocks)*100 
cat(pct_cpg, "blocks overlap with UCSC cpg Islands") 

### graph imputed methylation by sample####
pos <- names(index)[which(index$block == names(blocks[2]))]
boxplot(imp_meth[pos,], main="Sample variation block 2") 

