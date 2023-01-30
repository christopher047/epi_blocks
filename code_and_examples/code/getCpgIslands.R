
### cpg #####a
### sample overlap ###

getCpgIslands <- function(infile="../reference/cpgIslandExtV36.txt") 
	{
	cpg <- read.table(infile, sep="\t")
        colnames(cpg) <- c("bin", "chr", "start", "end", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
        cpg       <- GRanges(cpg)
        common    <- paste0("chr", c(1:22,"X", "Y"))
        cpg       <- keepSeqlevels(cpg, common, pruning.mode="coarse")
	return(cpg) 
	}


