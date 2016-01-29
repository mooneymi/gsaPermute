snpsToGenes = function(snpLocs, geneLocs, dist=0) {
	## snpLocs = 3 column dataframe: snpName, chr, position
	## geneLocs = 4 column dataframe: geneName, chr, start, end
	## dist = size of window around gene in which a SNP can be located (number of bases)
	snps2genes = list()
	for (i in 1:dim(geneLocs)[1]) {
		chr = geneLocs[i,2]
		start = as.numeric(geneLocs[i,3]) - dist
		end = as.numeric(geneLocs[i,4]) + dist
		snpList = snpLocs[snpLocs[,2] == chr & snpLocs[,3] >= start & snpLocs[,3] <= end, 1]
		if (length(snpList) > 0) {
			snps2genes[[ geneLocs[i,1] ]] = snpList
		}
	}
	return(snps2genes)
}


getGeneSizes = function(snpGeneMap) {
	## Count the number of SNPs assigned to each gene and create a named list containing the size of each gene
	snpCounts = lapply(snpGeneMap, function(x){length(x)})	
	return(snpCounts)
}


makeGeneBins = function(geneSnpCounts, minBinSize=25) {
	numGenes = length(geneSnpCounts)
	## Check that the number of available genes is not less than the bin size
	if (numGenes < 0.75*minBinSize) {
		e = simpleError('The minimum bin size is greater than the number of available genes.')
		stop(e)
	}	
	maxGeneSize = max(unlist(geneSnpCounts))
	
	## Create a named list containing vectors of genes grouped by size (the number of SNPs assigned to the gene)
	bins = as.list(rep(NA,maxGeneSize))
	names(bins) = c(1:maxGeneSize)
	for (i in 1:maxGeneSize) {
		bins[[i]] = names(geneSnpCounts[geneSnpCounts == i])
	}
	
	bounds = matrix(NA, nrow=0, ncol=2)
	i = 1
	while (i <= maxGeneSize) {
		start = i
		end = i
		count = length(bins[[start]])
		## Check the bin size and add genes to create a bin of at least the minimum size
		while (count < minBinSize & end <= (maxGeneSize-1)) {
			end = end + 1
			count = count + length(bins[[end]])
			
			## Check that there are enough genes remaining (75% of the minimum bin size), otherwise add all genes to current bin
			tmpIdx = end + 1
			tmpCount = length(bins[[tmpIdx]])
			while (tmpCount < (0.75*minBinSize)) {
				if (tmpIdx < maxGeneSize) {
					tmpIdx = tmpIdx + 1
					tmpCount = tmpCount + length(bins[[tmpIdx]])
				} else {
					count = count + tmpCount
					end = tmpIdx
					break
				}
			}
		}
		i = end + 1
		bounds = rbind(bounds, c(start, end))
	}
	
	return(list(simpleBins=bins, binBounds=bounds))
}


getGeneSetComposition = function(geneSnpCounts, geneSet) {
	snpCounts = c()
	for (gene in geneSet) {
		snpCounts = c(snpCounts, geneSnpCounts[[gene]])
	}
	return(table(snpCounts))
}


makeRandomGS = function(geneSetComp, geneBins) {
	bins = geneBins$simpleBins
	bounds = geneBins$binBounds
	
	randGS = c()
	sizes = names(geneSetComp)
	## Randomly select the correct number of genes of each size to match the composition of the specified gene set 
	for (size in sizes) {
		geneCount = geneSetComp[[size]]
		size = as.numeric(size)
		curBounds = bounds[bounds[,1] <= size & bounds[,2] >= size,]
		availGenes = as.character(unlist(bins[curBounds[1]:curBounds[2]]))
		availGenes = availGenes[!availGenes %in% randGS]
		randGS = c(randGS, sample(availGenes, geneCount))
	}
	return(randGS)
}


