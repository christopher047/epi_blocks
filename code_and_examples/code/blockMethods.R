
ReadRawFiles <- function(object, input_dir="./", achr="chr17")
                {
		#### reads raw coverage and raw methylation files ####	
		require(data.table) 	
		my_read <- function(filename) 
			{
			ff <- data.frame(fread(filename,sep="\t"))
                	rownames(ff) <- ff$V1
                	ff <- ff[,-1]
                	return(ff)
			}
		cat(paste("running", achr,"\n"))  
		object@achr	    <- achr
		object@raw_cov      <- my_read(paste0(input_dir, achr, "_tot_cov.tab"))  
        	object@raw_meth     <- my_read(paste0(input_dir, achr, "_tot_meth.tab")) 
		return(object) 
                }

createIndex  <- function(object) 
		{
		######## paramaters #####	
		tot_cov  <- object@raw_cov
		tot_meth <- object@raw_meth 
		#### keep coverage in index ####
		cat("creating index and row coverage\n") 
		chr            <- gsub("\\..*$","", rownames(tot_cov))
                pos            <- gsub("^.*\\.","", rownames(tot_cov))
                index          <- GRanges(data.frame("chr"=chr, "start"=pos, "end"=pos))
		names(index)   <- paste0(chr, ".", pos) 
        	########### create index ########
        	index            <- index[rownames(tot_cov)]
		index            <- sort(index)
		index$site_reads <- rowSums(tot_cov[names(index),], na.rm=T)  
        	index$site_mean  <- rowMeans(tot_cov[names(index),], na.rm=T) 	
		############ make sure they are in the same order ##########
        	object@raw_cov  <- tot_cov[names(index),]
        	object@raw_meth <- tot_meth[names(index),]
		object@index    <- index 
		return(object) 
		}

getSegs <- function(object, min.seg=5)
        {
	#### paramaters ####
	index     <- object@index
	#### get closest cutoff to 3000 ####
        test  <- diff(start(index))
        keep <- as.numeric(gsub("%", "", (names(which.max(which(quantile(test, seq(0,1,0.01)) < 3000))))))/100
        cutoff <- as.integer(quantile(test, keep))
        cat(paste("segment cutoff is", cutoff, "\n"))
        cat(paste("quantile is", keep, "\n"))
        #setup segmentation list
        test <- c(0, test)
        test[test >  cutoff] <- 0
        test <- setNames(test, as.character(names(index)))
        segs <- vector(mode="list", sum(test==0))
        #segment change to GRanges
        j = 0
        #index <- c()
        for(i in 1:length(test))
                {
                if (test[i] == 0){j <- j+1}
                segs[[j]] <- c(segs[[j]], names(test)[i])
                }
        #keep only groups of min.block or more
        segs <- lapply(segs, na.omit)
        segs <- segs[which(as.numeric(lapply(segs, length)) >=min.seg)]
	object@segs <- segs 
	return(object) 
        }


getGRangeIndex <- function(xx, achr)
        {
	### this takes chromosome, smallest, largest to create a GRange	
        st    <- lapply(xx, function(i){min(as.integer(gsub(".*\\.", "", i)))})
        ed    <- lapply(xx, function(i){max(as.integer(gsub(".*\\.", "", i)))})
        gg    <- GRanges(data.frame("chr"=achr, start=as.integer(st), end=as.integer(ed)))
        return(sort(gg))
        }

getSegRanges <- function(segs) 
	{	
	### this function converts segs to GRanges adds number of M sites as L	
	achr   <- gsub("\\..*", "", segs[[1]][1])
	temp   <- getGRangeIndex(segs, achr) 
        temp$L <- lapply(segs, length)
        return(temp)
	}


createSegIndex <- function(segs, achr) 
	{
	### this names range by chromosome and start position	
	names(segs)   <- paste0(achr, ".", 1:length(segs))
        sindex        <- getGRangeIndex(segs, achr)
        names(sindex) <- paste0(achr, ".", 1:length(segs))
	return(sindex) 
	}	

imputeKNNbySegs <- function(object, ncores=16, nb=5) 
	{
	cat("imputing by segment \n") 	
	require(impute)
	require(parallel) 
	segs <- object@segs 	
	temp <- object@raw_meth/object@raw_cov
        temp[temp == Inf] <- NA
        x1       <- mclapply(1:length(segs), function(i){impute.knn(t(temp[segs[[i]],]), k=nb)$data}, mc.cores=ncores)
        b0       <- data.frame(t(do.call("cbind", x1)))
	object@imp_mat <- b0
	return(object) 
	}

methylDist <- function(X,pos,bwd = 300)
        {
        d1 = cor(t(X))
        d1[is.na(d1)] <- 1
        d2 <- exp(-as.matrix(dist(pos)^2)/(2*bwd^2))
        d <- 1-(d1*d2)
        return(d)
        }

endreCluster <- function(bb=bb, bwd=300, min.block=3, hclust=0.4)
        {
        pos    <- as.integer(unlist(lapply(strsplit(rownames(bb),split="\\."), "[[", 2)))
        d      <- methylDist(bb,pos, bwd = bwd)
        clRes  <- hclust(as.dist(d),method = 'single')
        cutRes <- cutree(clRes,h = hclust)
        tab    <- table(cutRes)
        temp   <- aggregate(names(cutRes), by=list("members"=as.integer(cutRes)), c)
        temp   <- temp[temp$member %in% names(which(tab >=min.block)),]$x
        if(is.matrix(temp)) {temp <- list(as.vector(temp))}
        if(length(temp) == 0){return(c())}
        return(temp)
        }

getBlocks <- function(object, bwd=300, min.block=5, hclust=0.4, ncores=16)
        {
	#### paramaters ####
	b0    <- object@imp_mat
        achr  <- object@achr
	segs  <- object@segs
	index <- object@index	
        #### find blocks from segments ####     
        t0        <- mclapply(1:length(segs), function(i){endreCluster(bb=b0[segs[[i]],], bwd=bwd, min.block=min.block, hclust=hclust)}, mc.cores=ncores)
        names(t0) <- names(segs)
        t0        <- unlist(t0, recursive=F)
        t0        <- t0[which(as.numeric(lapply(t0, length)) > 0)]
        if (length(t0) < 1)
		{
		object@blocks <- GRanges() 
		object@index  <- GRanges() 
		object@imp_mat <- data.frame()  
		return(object) 
		}		
	### what if there are no blocks ###
        ### create block GRanges and update index###
        blocks         <- getGRangeIndex(t0, achr)
        names(blocks)  <- paste0(achr,".", 1:length(blocks))
        ov             <- findOverlaps(index, blocks)
        index$block    <- NA
        index[queryHits(ov)]$block <- names(blocks[subjectHits(ov)])
        index <- index[!is.na(index$block),]
	#### return ####
	object@blocks <- blocks 
	object@index  <- index 
	return(object) 
        }

filterIQR <- function(object, iqr_cutoff=10)
        {
	if (length(object@blocks) < 1){return(object)}
	#### paramaters ####
	b0     <- object@imp_mat 
	index  <- object@index 
	blocks <- object@blocks 
	##### remove outliers ####	
        temp_imp      <- b0[names(index),]
        xx            <- split(temp_imp, index$block)
        xx            <- lapply(xx, as.matrix)
        iqr           <- lapply(xx, IQR)
        rr            <- lapply(lapply(xx, range), diff)
        fx            <- as.numeric(rr)/as.numeric(iqr[names(rr)])
        fx            <- setNames(fx, names(xx))
        fx            <- fx[is.finite(fx)]
        fx            <- fx[which(fx <= iqr_cutoff)]
	keep          <- names(fx) 
	#### check to see if we have no blocks 
	if (length(keep) < 1) 
		{
		object@blocks <- GRanges() 
		object@index <- GRanges() 
		}
	if (length(keep)>0)
		{	
		object@blocks <- sort(blocks[keep])
        	object@index  <- index[index$block %in% keep]
		}
	return(object) 
        }

calcBlockCoverage <- function(object, min.block=5, ncores=16) 
	{
	#### wei error ###
	if (length(object@blocks) < 1){return(object)}
	#### paramaters ####
	blocks          <- object@blocks
        index           <- object@index
	b0              <- object@imp_mat 
	cov             <- object@raw_cov
	b0              <- b0[names(index),]
	cov             <- cov[names(index),]
	cov[is.na(cov)] <- 0 
	#### initialize ####
 	blocks$block_mean <- 0
	blocks$block_sum  <- 0
	blocks$L          <- 0
	blocks$PC	  <- 0 
	blocks$imp        <- 0 
	#### calculate average methylation ###
	xx <- split(b0, index$block) 
	PC <- lapply(lapply(xx, colMeans), mean)
	blocks[names(PC)]$PC <- as.numeric(PC)	
	#### calculate number of cpgs  ####
	xx <- split(cov, index$block) 
	L  <- lapply(xx, nrow) 
	blocks[names(L)]$L <- as.numeric(L) 
	#### calculate total reads per block ####
	block_sum <- lapply(xx, sum)
	blocks[names(block_sum)]$block_sum <- as.numeric(block_sum)  
	#### calculate average reads per block ####
	block_mean <- lapply(lapply(xx, rowSums), mean) 
	blocks[names(block_mean)]$block_mean <- as.numeric(block_mean) 
	#### calculate impution ####
	imp <- lapply(xx, function(xx){sum(xx==0)/(nrow(xx)*ncol(xx))})	
	blocks[names(imp)]$imp <- as.numeric(imp) 
	### verify and update####
	blocks <- blocks[which(blocks$L >= min.block)]
        index  <- index[index$block %in% names(blocks)]
        ### sort and update opject 
	index <- sort(index) 
	object@imp_mat <- b0[names(index),]
	object@raw_cov <- cov[names(index),]
	object@blocks  <- blocks
        object@index   <- index
        return(object)
	}

getReturnData <- function(object)
        {
	### what if there are no blocks 	
        b0               <- object@imp_mat
	blocks           <- object@blocks 
	index            <- object@index
	raw_cov          <- object@raw_cov
	segs             <- object@segs
	cat(paste("found", length(blocks), "blocks\n"))  
	return(list(imp_meth=b0, blocks=blocks, index=index, raw_cov=raw_cov, segs=segs)) 	
        }



