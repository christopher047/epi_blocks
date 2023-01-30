library(GenomicRanges)
source("blockMethods.R")

### class definition  ###

setClass("meth", slots=list(achr    = "character",  #the chromosome being processed
                            imp_mat ="data.frame",  #the imputed matrix 
                            raw_cov = "data.frame", #the raw coverage 
                            raw_meth = "data.frame",#the raw methylation
                            segs    ="list",        #the segments 
                            index   ="GRanges",     #GRange index of each site
                            blocks  ="GRanges",     #GRange index of each block
                            readFiles ="function",  #reads the files 
                            createIndex = "function", #create index  
                            getSegs   = "function", #divided chromosomes into segments 
                            imputeKNNbySegs = "function",#imputes the imp_matrix by segment, imp_mat is raw_meth/raw_counts, raw_count =0 gives NA 
                            getBlocks       = "function", #get blocks by variation and distance
                            filterIQR       = "function", #filters blocks with one or two outliers 
                            calcBlockImputation = "function", #calculates the percent impution per block
                            calcBlockCoverage   = "function", #calculates L(number of cpgs), block_sum(number total reads), block_mean(mean relative methylation)
                            getReturnData       = "function"))# returns blocks and imp_mat of blocks 


### initialize class  ####

createMethObject <- function()
        {
        temp <- new("meth")
        temp@readFiles           <- ReadRawFiles
        temp@createIndex         <- createIndex
        temp@getSegs             <- getSegs
        temp@imputeKNNbySegs     <- imputeKNNbySegs
        temp@getBlocks           <- getBlocks
        temp@filterIQR           <- filterIQR
        temp@calcBlockCoverage   <- calcBlockCoverage
        temp@getReturnData       <- getReturnData
        return(temp)
        }

### intialize instance  ###

runChromosome <- function(achr, input_dir="../chrs/", min.seg=20, bwd=300, min.block=10, hclust=0.4, iqr_cutoff=10, nb=5, ncores=16)
        {
        test <- createMethObject()
        test <- test@readFiles(test, input_dir=input_dir, achr=achr)
        test <- test@createIndex(test)
        test <- test@getSegs(test, min.seg=min.seg)
        test <- test@imputeKNNbySegs(test, ncores=ncores, nb=nb)
        test <- test@getBlocks(test, bwd=bwd, min.block=min.block, hclust=hclust, ncores=ncores)
        test <- test@filterIQR(test, iqr_cutoff = iqr_cutoff)
        test <- test@calcBlockCoverage(test, min.block=min.block, ncores=ncores)
        res  <- test@getReturnData(test)
        if(!all(!is.na(test@index))){break}
        return(res)
        }

