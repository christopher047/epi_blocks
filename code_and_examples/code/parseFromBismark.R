
getSafeName <- function(infile) 
	{
	temp <- data.frame(matrix(nrow=1, ncol=1, 0, dimnames=list("a", basename(infile))))
	return(colnames(temp)[1]) 
	}

getCovMeth   <- function(infile, n=4, output_dir="./test/")
	{
	require(data.table) 
	require(tools) 
	my_start <- Sys.time() 
	cat(paste("processing", infile, "\n")) 
	out_cov      <- paste0(output_dir, "/", "COV_" , getSafeName(infile))   
	out_meth     <- paste0(output_dir, "/", "METH_", getSafeName(infile))
	#if (file_ext(infile) == "gz"){infile <- paste("zcat", infile)}
	df           <- fread(infile, sep="\t")
	colnames(df) <- c("chr", "start", "end", "meth_pct", "methylated", "not_methylated")
	df           <- df[df$methylated + df$not_methylated > n]
	df           <- df[grep("chr", df$chr)]
	df[, id := do.call(paste, c(.SD, sep=".")), .SDcols=c(1,2)]
	fwrite(data.table(df$id, df$methylated +df$not_methylated), out_cov, sep="\t")  
	fwrite(data.table(df$id, df$methylated), out_meth, sep="\t")
	cat(paste("finished", Sys.time() - my_start, "\n")) 
	return(df$id) 
	}

getInitialMatrix <- function(files, n=4, n_imp=70, output_dir="./test/")
       	{
	cat("creating master matrix\n")
	all_ids        <- lapply(files, getCovMeth, n=n, output_dir=output_dir) 
	files          <- basename(list.files(path=output_dir, pattern="^COV_", full.name=T)) 
	files          <- gsub("COV_", "", files)
	names(all_ids) <- files 
	uu             <- Reduce(union, all_ids) 
	### create matrix for all samples ####
	dummy          <- data.frame(matrix(nrow=length(uu), ncol=length(all_ids), 0, dimnames=list(uu, names(all_ids))))
        for(i in 1:length(all_ids)) { dummy[all_ids[i][[1]], names(all_ids)[i]] <- 1}
	### make a graph to show users samples with low values ####
	samp_id <- data.frame(colSums(dummy)/nrow(dummy)*100, colSums(dummy))
	colnames(samp_id) <- c(paste("pc of total sites ncov gt", n), paste("total sites gt", n)) 
	write.table(samp_id, paste0(output_dir, "sample_overview.tab"), sep="\t") 
	### remove sites less than n_imp percentage ####
	dd   <- rowSums(dummy)/ncol(dummy)*100
	dummy <- dummy[which(dd >= n_imp),]
        dummy[dummy==1] <- 0
	cat("finished\n") 
	return(dummy)
	}	

createCoverageMatrix <- function(dummy_mat, output_dir="./test/", pattern="^COV_") 
	{
	cat(paste("creating", pattern,  "total matrix\n"))  	
	files <- list.files(path=output_dir, pattern=pattern, full.name=T)
	stub  <- gsub(pattern, "", basename(files))
	for(i in 1:length(files)) 
		{
		temp  <- fread(files[i], sep="\t") 
		temp <- temp[temp$V1 %in% rownames(dummy_mat),]
		dummy_mat[temp$V1, stub[i]] <- as.integer(temp$V2) 
		}
	cat("finished\n") 
	return(dummy_mat) 
	}

splitIntoChromosomes <- function(count_mat, output_dir="./test/", prefix="_tot_cov.tab")
	{	
	require(GenomicRanges) 	
	cat("splitting into chromosomes\n") 	
	chr     <- gsub("\\..*", "", rownames(count_mat))
	pos     <- gsub(".*\\.", "", rownames(count_mat)) 
	temp    <- GRanges(data.table(chr=chr, start=pos, end=pos, id=row.names(count_mat)))
	temp    <- sort(temp)
	all_chr <- unique(as.character(seqnames(temp)))
	for (i in 1:length(all_chr)) 
		{
		output_file <- paste0(output_dir, all_chr[i],prefix)  	
		write.table(count_mat[temp[which(seqnames(temp) == all_chr[i])]$id,], output_file, sep="\t") 	
		}	
	}

clean_up <- function(output_dir="./test/") 
	{
	cov_files <- list.files(path=output_dir, pattern="^COV_", full.name=T) 
	meth_files <- list.files(path=output_dir, pattern="^METH_", full.name=T)
	sapply(cov_files, unlink)
	sapply(meth_files, unlink)
	}


parseBismark <- function(files, n=4, n_imp=70, output_dir="../chrs/") 
	{
	my_start            <- Sys.time()
	dir.create(output_dir, showWarnings = F) 
	dummy               <- getInitialMatrix(files, n=n, n_imp=n_imp, output_dir=output_dir) 
	tot_cov             <- createCoverageMatrix(dummy, output_dir=output_dir, pattern="^COV_") 
	tot_meth            <- createCoverageMatrix(dummy, output_dir=output_dir, pattern="^METH_") 
	tot_cov[tot_cov==0] <- NA
	splitIntoChromosomes(tot_cov, output_dir=output_dir, prefix="_tot_cov.tab") 
	splitIntoChromosomes(tot_meth, output_dir=output_dir, prefix="_tot_meth.tab") 
	clean_up(output_dir=output_dir) 
	cat(paste("\n\nTOTAL", Sys.time() - my_start, "\n"))
	}

