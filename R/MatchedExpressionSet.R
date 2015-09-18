# TODO: Add comment
# 
# Author: jthoeve
###############################################################################

matchEsets <- function(...) {
	new('MatchedExpressionSet', ...)
}

setMethod('initialize', 'MatchedExpressionSet', function(.Object, ...) {
	
	# Do argument checks!
	args <- list(...)

	if ( length(args) <= 1 ) {
		stop("Provide a least 2 *uniquely* named ExpressionSets!")
	}
	
	if ( !all(sapply(args, inherits, what='ExpressionSet')) ) {
		stop("Provide a least 2 *uniquely* named ExpressionSets!")
	}
	
	if ( length(unique(names(args))) != length(names(args))) {
		stop("Provide a least 2 *uniquely* named ExpressionSets!")
	}
	
	if ( !all(nchar(names(args)) >= 1 ) ) {
		stop("Provide a least 2 *uniquely* named ExpressionSets!")
	}
	matched.samples <- sampleNames(args[[1]])
	for (i in 2:length(args)) {
		matched.samples <- intersect(matched.samples, sampleNames(args[[i]]))
	}
	
	.Object@expressionSets <- lapply(args, function(eset) eset[, matched.samples])
	message('Matched ', length(matched.samples), ' samples and combined ', length(args), ' ExpressionSets (', paste(names(args), collapse=', '),').')


	# Calculating correlation matrix
	vals1 <- exprs(.Object@expressionSets[[1]])
	vals2 <- exprs(.Object@expressionSets[[2]])

	message('Calculating correlation matrix [', nrow(vals1), ', ', nrow(vals2), ']...')
	.Object@correlationMatrix <- cor(t(vals1), t(vals2), method='pearson', use='pairwise.complete.obs')
	names(dimnames(.Object@correlationMatrix)) <- names(.Object@expressionSets)


	return(.Object)
})

setMethod('show', 'MatchedExpressionSet', function(object) {
	cat('An object of class "MatchedExpressionSet"\n')
	cat('ExpressionSets:', names(object@expressionSets), '\n')
	show(object@expressionSets)

	cat('Phenodata:\n')
	str(pData(object))

	cat("\nInteraction DataFrames (", length(object@interactionDataFrames), '):\n',
	paste('  ', names(object@interactionDataFrames), ' (', sapply(object@interactionDataFrames, nrow), ')', sep='', collapse='\n'), '\n', sep='')

	cat("\nThresholds:", thresholds(object), '\n')
	cat("\nName:", name(object), '\n')	
})

setMethod('name', 'MatchedExpressionSet', function(object) {
	object@name
})

setMethod('thresholds', 'MatchedExpressionSet', function(object) {
  object@thresholds
})

setReplaceMethod('name', 'MatchedExpressionSet', function(object, value) {
	object@name <- value
	object
})

setMethod('sampleNames', 'MatchedExpressionSet', function(object) {
	sampleNames(object@expressionSets[[1]])
})

setMethod('pData', 'MatchedExpressionSet', function(object) {
	pData(object@expressionSets[[1]])
})

setMethod('eset', c('MatchedExpressionSet', 'character'), function(object, id) {
	object@expressionSets[[id]]
})

setMethod('esetNames', c('MatchedExpressionSet'), function(object) {
	names(object@expressionSets)
})


setMethod('subset', c('MatchedExpressionSet'), function(object, i, j, thresholds=NULL, name, removeInteractions=FALSE) {

#	stop('Indexing for matchedEsets is not implemented!')

	if ( missing(i) )
		i <- NULL

	if ( missing(j) )
		j <- NULL

	if (!is.null(thresholds))
    	stopifnot(length(thresholds) == length(object@expressionSets))

	stopifnot(!missing(name))

	subsettedexpressionSets <- lapply(names(object@expressionSets), function(eset.id) {

		eset <- object@expressionSets[[eset.id]]

		if ( is.null(i) ) {
			idx.i <- featureNames(eset)
		} else if ( is.null( i[[eset.id]] ) ) {
			idx.i <- featureNames(eset)
		} else {
			i <- i[[eset.id]]
			idx.i <- i[ i %in% featureNames(eset) ]
		}

		if ( is.null(j) ) {
			idx.j <- sampleNames(eset)
		} else {
			idx.j <- j[ j %in% sampleNames(eset) ]
		}


		sub.eset <- eset[idx.i, idx.j]

		if (!is.null(thresholds[[eset.id]])) {
		  means <- apply(exprs(sub.eset)[, , drop=FALSE], 1, mean, na.rm=TRUE)
		  idx   <- means >= thresholds[[eset.id]]
		  idx[is.na(idx)] <- FALSE
      
	      message('Filtering on ', eset.id, ', mean >= ', thresholds[[eset.id]], ': ', sum(idx), '/', length(idx))
		  sub.eset <- sub.eset[idx, , drop=FALSE]
		}
		
		sub.eset
    
 	})

	names(subsettedexpressionSets) <- names(object@expressionSets)

	fun <- function(...) {
		new('MatchedExpressionSet', ...)
	}

	newMes <- do.call('fun', subsettedexpressionSets)

	if (!removeInteractions) {
		for (id in interactionNames(object)) {
			interactions(newMes, id) <- interactions(object, id)
		}
	}
	name(newMes) <- name
	
	if (!is.null(thresholds)) {
		newMes@thresholds <- thresholds
	}

	newMes
})

setMethod('subsetByFeatures', c('MatchedExpressionSet', 'list'), function(object, indexList) {

	# expressionSets
	for (id in names(indexList)) {
		object@expressionSets[[id]] <- object@expressionSets[[id]][indexList[[id]], ]
	}

	# interactionDataFrames = 'list',
	for(name in interactionNames(object)) {
		idf <- object@interactionDataFrames[[name]]

		idx1 <- idf[[1]] %in% featureNames(object@expressionSets[[1]])
		idx2 <- idf[[2]] %in% featureNames(object@expressionSets[[2]])

		object@interactionDataFrames[[name]]  <- idf[idx1 & idx2, ]
	}

	# correlationMatrix
	idx1 <- rownames(object@correlationMatrix) %in% featureNames(object@expressionSets[[1]])
	idx2 <- colnames(object@correlationMatrix) %in% featureNames(object@expressionSets[[2]])

	object@correlationMatrix <- object@correlationMatrix[idx1, idx2, drop=FALSE]
	object
})



setMethod('[', c('MatchedExpressionSet', 'character', 'character'), function(x, i, j, drop=FALSE) {


	stop('Indexing for matchedEsets is not implemented!')

	if ( missing(i) )
		i <- NULL

	if ( missing(j) )
		j <- NULL


	subsettedexpressionSets <- lapply(names(x@expressionSets), function(eset.id) {

		eset <- x@expressionSets[[eset.id]]

		if ( is.null(i) ) {
			idx.i <- featureNames(eset)
		} else if ( is.null( i[[eset.id]] ) ) {
			idx.i <- featureNames(eset)
		} else {
			i <- i[[eset.id]]
			idx.i <- i[ i %in% featureNames(eset) ]
		}

		if ( is.null(j) ) {
			idx.j <- sampleNames(eset)
		} else {
			idx.j <- j[ j %in% sampleNames(eset) ]
		}

		eset[idx.i, idx.j]
	})
	names(subsettedexpressionSets) <- names(x@expressionSets)

	fun <- function(...) {
		new('MatchedExpressionSet', ...)
	}

	newMes <- do.call('fun', subsettedexpressionSets)

	newMes@interactionDataFrames <- x@interactionDataFrames

	newMes
})


setMethod('interactionNames', c('MatchedExpressionSet'), function(object) {
	names(object@interactionDataFrames)
})

setMethod('interactions', c('MatchedExpressionSet', 'character'), function(object, id) {
	object@interactionDataFrames[[id]]
})

setReplaceMethod('interactions', 'MatchedExpressionSet', function(object, id, value, showWarnings=TRUE) {

	if (is.null(value)) {
		object@interactionDataFrames[[id]] <- NULL
	} else {

		if (!all(value[[1]] %in% featureNames(object@expressionSets[[1]]))) {
			if (showWarnings) warning("Rows in 'values' where column 1 does not match featureNames of eset 1 are left out.")
		}
		idx1 <- value[[1]] %in% featureNames(object@expressionSets[[1]])
		

		if (!all(value[[2]] %in% featureNames(object@expressionSets[[2]]))) {
			if (showWarnings) warning("Rows in 'values' where column 2 does not match featureNames of eset 2 are left out.")
		}
		idx2 <- value[[2]] %in% featureNames(object@expressionSets[[2]])
		
		object@interactionDataFrames[[id]] <- value[idx1 & idx2, ]
	}

	object
})

setMethod('interactionCounts', c('MatchedExpressionSet', 'list'), function(object, genesets, split=FALSE, return.genes=FALSE, negative.correlations.only=FALSE) {
	
	ids <- interactionNames(object)
	ics <- lapply(ids, function(id) {
		idf <- interactions(object, id)

		if (negative.correlations.only) {
			idf$cors <- correlations(object, id)
			idf <- idf[idf$cors < 0, ]
		}

		sapply(genesets, function(genes) {
			if (split) {
				splitlist <- split(idf[[2]], idf[[1]])
				sapply(splitlist, function(genes2) {
					if (return.genes) {
						paste(intersect(genes, genes2), collapse=', ')
					} else {
						sum(genes %in% genes2)
					}
				})
			} else {
				if (return.genes) {
					paste(intersect(genes %in% idf[[2]]), collapse=', ')
				} else {
					sum(genes %in% idf[[2]])
				}
			}
		})
	})
	names(ics) <- ids
	ics
})








setMethod('correlationMatrix', 'MatchedExpressionSet', function(object) {
	object@correlationMatrix
})

setMethod('correlations', 'MatchedExpressionSet', function(object, id, na.rm=FALSE, invert=FALSE) {

	idxOfPairsInMatrix <- function(x, y, mat) {

		if (missing(mat)) {
			ux <- sort(unique(x))
			uy <- sort(unique(y))

			mat <- matrix(FALSE, nrow=length(ux), ncol=length(uy))
			rownames(mat) <- ux
			colnames(mat) <- uy
		}

		ix <- match(x, rownames(mat))
		iy <- match(y, colnames(mat))

		nrow(mat) * (iy-1) + ix 
	}

	cm <- object@correlationMatrix

	if (missing(id)) {
		cors <- as.vector(cm)
	} else {
		idf <- interactions(object, id)


		idx <- idxOfPairsInMatrix(idf[[1]], idf[[2]], cm)
		if(invert) {
			cors <- cm[-idx]
		} else {
			cors <- cm[idx]
		}
	}

	if (na.rm) {
		cors <- cors[!is.na(cors)]
	}

	return(cors)
})









# mapTo : f.e.'mirbase'
# mapping : dataframe containing two columns
#
setMethod('mapEset', c('MatchedExpressionSet'), function(object, eset.id, mapTo, mapping, name, method=c('mean', 'maxvar')) {

	mapFrom <- annotation(eset(mEset, eset.id))
	stopifnot(mapFrom != mapTo)

	method <- match.arg(method)

	eset <- object@expressionSets[[eset.id]]


	#uniqueTos <- unique(mapping[[mapTo]])

	mapping <- mapping[mapping[[mapFrom]] %in% featureNames(eset), ]
	
	mapList <- split(mapping[[mapFrom]], mapping[[mapTo]])


	mappedExprs <- t(sapply(names(mapList), function(toID) {
		fromIDs <- mapList[[toID]]

		idx <- match(featureNames(eset), fromIDs)
		idx <- idx[!is.na(idx)]

		stopifnot(length(idx) != 0)

		if (length(idx) == 1) {
			exprs(eset)[fromIDs[idx], ]
		} else {

			if (method == 'maxvar') {
				vars <- apply(exprs(eset)[fromIDs[idx], , drop=FALSE], 1, var)
				exprs(eset)[fromIDs[idx[which.max(vars)]], ]
			} else if (method == 'mean')
				colMeans(exprs(eset)[fromIDs[idx], ], na.rm=TRUE)
			else {
				stop(paste('method: "', method, '" not implemented'))
			}
		}
	}))

	expressionSets <- object@expressionSets
	exprs(expressionSets[[eset.id]]) <- mappedExprs
	annotation(expressionSets[[eset.id]]) <- mapTo

	fun <- function(...) {
		new('MatchedExpressionSet', ...)
	}

	newMes <- do.call('fun', expressionSets)

	newMes@thresholds <- object@thresholds
	newMes@name <- object@name

	newMes
})








setMethod('abundancePlot', 'MatchedExpressionSet', function(
	object,
	id,
	fixed.eset               = c('mirna', 'mrna'),
	fixed.assayDataElement   = 'exprs',
	fixed.threshold          = NA,
	varying.assayDataElement = 'exprs',   # 'A'
	n.varying.thresholds     = 20,
	y = c('shift', 'p-value', 'n.varying'),
	...) {




	fixed.eset   <- match.arg(fixed.eset)
	y            <- match.arg(y, several.ok=TRUE)

	# opposite of the fixed.eset
	varying.eset <- switch(fixed.eset,
		mirna = 'mrna',
		mrna  = 'mirna'
	)


	# 
	fixed.means <- apply(assayDataElement(eset(mEset, fixed.eset), fixed.assayDataElement), 1, mean, na.rm=TRUE)

	if (!is.na(fixed.threshold)) {
		fixed.idx <- fixed.means >= fixed.threshold
		fixed.idx[is.na(fixed.idx)] <- FALSE
	} else {
		fixed.idx <- rep(TRUE, length(fixed.means))
	}

	
	varying.means <- apply(assayDataElement(eset(mEset, varying.eset), varying.assayDataElement), 1, mean, na.rm=TRUE)

	if (varying.eset == 'mirna' & varying.assayDataElement == 'exprs') {
		varying.thresholds <- 10^(seq(-7, -1, .25))

	} else {

		varying.thresholds <- seq(min(varying.means, na.rm=TRUE), max(varying.means, na.rm=TRUE), length.out=n.varying.thresholds)
		varying.thresholds <- varying.thresholds[1:(n.varying.thresholds-1)]
	}

	### Loop over thresholds

	result <- sapply(varying.thresholds, function(threshold) {

		varying.idx <- varying.means >= threshold
		varying.idx[is.na(varying.idx)] <- FALSE


		message('threshold ', threshold, ': ', sum(varying.idx), ' ', varying.eset, "'s...")


		subsetList <- list(fixed.idx, varying.idx)
		names(subsetList) <- c(fixed.eset, varying.eset)

		subsettedMatchedEset <- subsetByFeatures(mEset, subsetList)


		cors.bg <- correlations(subsettedMatchedEset, id, na.rm=TRUE, invert=TRUE) 
		cors.fg <- correlations(subsettedMatchedEset, id, na.rm=TRUE, invert=FALSE) 

		p.val <- wilcox.test(cors.fg, cors.bg, alternative='less')$p.value


		c(
			shift          = median(cors.bg) - median(cors.fg),
			'p-value'      = -log10(p.val),
			n.varying      = sum(varying.idx)
		)
	})


	old.par <- par(mfrow=n2mfrow(length(y)))
	for (var in y) {

		if (varying.eset == 'mirna' & varying.assayDataElement == 'exprs') {
			plot(varying.thresholds, result[var, ], log='x', xaxt='n', type='l', main=paste('Abundance plot - fixed ', fixed.eset, ' (', fixed.assayDataElement, ') at ', fixed.threshold, sep=''), xlab=paste('Threshold on ', varying.eset, ' (', varying.assayDataElement, ')', sep=''), ylab=var)
			axis(1, c(.1, .01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001))
		} else {
			plot(varying.thresholds, result[var, ], type='l', main=paste('Abundance plot - fixed ', fixed.eset, ' (', fixed.assayDataElement, ') at ', fixed.threshold, sep=''), xlab=paste('Threshold on ', varying.eset, ' (', varying.assayDataElement, ')', sep=''), ylab=var)
		}
	}
	par(old.par)

})









setMethod('kmPlot', 'ExpressionSet', function(object, feature, time='', status='', weights=1) {


	stopifnot(length(feature) == length(weights))
	stopifnot(all(feature %in% rownames(exprs(object))))


	X <- colSums(exprs(object)[feature, , drop=FALSE] * weights, na.rm=TRUE)
	group <- as.factor(as.numeric(X >= median(X)) + 1)
	levels(group) <- c('< median', '>= median')




	# COX model
	# Assume right sensored survival data

	Y <- pData(object)[, c(time, status)]

	na.idx <- apply(Y, 1, function(x) any(is.na(x)))

	Y <- Y[!na.idx, ]
	X <- X[!na.idx]
	group <- group[!na.idx]

	df <- data.frame(
		time   = Y[[time]],
		status = Y[[status]],
		x = group
	)

	feature.surv <- survfit(Surv(time, status) ~ x, data = df)
	plot(feature.surv, lty = 2:3, xlab=time, ylab=status)
	legend('bottomleft', levels(group), lty = 2:3) 
	title(paste("Kaplan-Meier Curves\nfor ", paste(feature, collapse=', '), sep='')) 
	lsurv2 <- survfit(Surv(time, status) ~ x, df, type='fleming') 
	plot(lsurv2, lty=2:3, fun="cumhaz", xlab=time, ylab="Cumulative Hazard")

	barplot(X, col=Y[, status])

})


setMethod('kmPlotPheno', 'MatchedExpressionSet', function(object, pheno, time='', status='') {

#	X <- rbind(
#		exprs(eset(object, id = 'mirna')),
#		exprs(eset(object, id = 'mrna'))
#	)[feature, ]
	

	group <- as.factor(pData(object)[, pheno])
	



	# COX model
	# Assume right sensored survival data

	Y <- pData(object)[, c(time, status)]

	na.idx <- apply(Y, 1, function(x) any(is.na(x))) | is.na(group)

	Y <- Y[!na.idx, ]
#	X <- X[!na.idx]
	group <- group[!na.idx]

	df <- data.frame(
		time   = Y[[time]],
		status = Y[[status]],
		x = group
	)

	pheno.surv <- survfit(Surv(time, status) ~ x, data = df)
	plot(pheno.surv, lty = 2:3, xlab=time, ylab=status)
	legend('bottomleft', levels(group), lty = 2:3) 
	title(paste("Kaplan-Meier Curves\nfor ", pheno, sep='')) 
	lsurv2 <- survfit(Surv(time, status) ~ x, df, type='fleming') 
	plot(lsurv2, lty=2:3, fun="cumhaz", xlab=time, ylab="Cumulative Hazard")

#	barplot(X, col=Y[, status])
})



setMethod('kmPlot', 'MatchedExpressionSet', function(object, feature, time='', status='') {

	X <- rbind(
		exprs(eset(object, id = 'mirna')),
		exprs(eset(object, id = 'mrna'))
	)

	if (!feature %in% rownames(X)) {
		warning(feature, ' not found in the data')
	} else {

		X <- X[feature, ]
		group <- as.factor(as.numeric(X >= median(X)) + 1)
		levels(group) <- c('< median', '>= median')




		# COX model
		# Assume right sensored survival data

		Y <- pData(object)[, c(time, status)]

		na.idx <- apply(Y, 1, function(x) any(is.na(x)))

		Y <- Y[!na.idx, ]
		X <- X[!na.idx]
		group <- group[!na.idx]

		df <- data.frame(
			time   = Y[[time]],
			status = Y[[status]],
			x = group
		)

		feature.surv <- survfit(Surv(time, status) ~ x, data = df)
		plot(feature.surv, lty = 2:3, xlab=time, ylab=status)
		legend('bottomleft', levels(group), lty = 2:3) 
		title(paste("Kaplan-Meier Curves\nfor ", feature, sep='')) 
		lsurv2 <- survfit(Surv(time, status) ~ x, df, type='fleming') 
		plot(lsurv2, lty=2:3, fun="cumhaz", xlab=time, ylab="Cumulative Hazard")

		barplot(X, col=Y[, status])
	}
})



setMethod('correlationBoxplot', 'MatchedExpressionSet', function(object, id, mir, genes, group='molecular_subtype', plot=TRUE) {
	mirExprs  <- exprs(eset(object, 'mirna'))
	mrnaExprs <- exprs(eset(object, 'mrna'))

	mrnaMeans <- apply(mrnaExprs, 1, var, na.rm=TRUE)

	if (is.null(genes)) {
		warning(mir, ' has no targets!')
	} else {

		cors.per.group <- lapply(split(sampleNames(object), pData(object)[[group]]), function(samples) {
			sapply(genes, function(gene) {
				cor(mirExprs[mir, samples], mrnaExprs[gene, samples], method='pearson', use='pairwise.complete.obs')
			})
		})
		if (plot) {
			boxplot(cors.per.group, main=paste(mir, '\n', paste(genes, collapse=', '), sep=''))
		}

		return(cors.per.group)
	}
})

if (FALSE) {
	setMethod('correlationBoxplot', 'MatchedExpressionSet', function(object, id, mir, cor.threshold, alternative=c('less', 'greater'), group='molecular_subtype') {
		
		alternative <- match.arg(alternative)

		mirExprs  <- exprs(eset(object, 'mirna'))
		mrnaExprs <- exprs(eset(object, 'mrna'))

		mrnaMeans <- apply(mrnaExprs, 1, var, na.rm=TRUE)

		mir <- match.arg(mir, rownames(mirExprs))

		idf <- interactions(object, id)
		idf$cors <- correlations(object, id)
		idf$mrna.mean <- mrnaMeans[idf[[2]]]
	#	idf <- idf[order(idf$mrna.mean, decreasing=TRUE), ]
		idf <- idf[order(idf$cor), ]
		idf <- switch(alternative,
			less    = idf[idf$cor < cor.threshold, ],
			greater = idf[idf$cor > cor.threshold, ]
		)


		idf.split <- split(idf[[2]], idf[[1]])
		genes <- idf.split[[mir]]
		if (is.null(genes)) {
			warning(mir, ' has no targets!')
		} else {

			#genes <- genes[1:min(top.n, length(genes))]
			cors.per.group <- lapply(split(sampleNames(object), pData(object)[[group]]), function(samples) {
				sapply(genes, function(gene) {
					cor(mirExprs[mir, samples], mrnaExprs[gene, samples], method='pearson', use='pairwise.complete.obs')
				})
			})

	#		boxplot(cors.per.group, main=paste(mir, '\n', paste(genes, collapse=', '), sep=''))
			boxplot(cors.per.group, main=paste(mir, '\n', length(genes), ' targets have correlation ', alternative,' than ', cor.threshold, sep=''))
		}
	})


	setMethod('correlationBoxplot', 'MatchedExpressionSet', function(object, id, mir, top.n=10, group='molecular_subtype') {
		mirExprs  <- exprs(eset(object, 'mirna'))
		mrnaExprs <- exprs(eset(object, 'mrna'))

		mrnaMeans <- apply(mrnaExprs, 1, var, na.rm=TRUE)

		mir <- match.arg(mir, rownames(mirExprs))

		# get the top.n genes.
		idf <- interactions(object, id)
		idf$cors <- correlations(object, id)
		idf$mrna.mean <- mrnaMeans[idf[[2]]]
	#	idf <- idf[order(idf$mrna.mean, decreasing=TRUE), ]
		idf <- idf[order(idf$cor), ]


		idf.split <- split(idf[[2]], idf[[1]])
		genes <- idf.split[[mir]]

		if (is.null(genes)) {
			warning(mir, ' has no targets!')
		} else {

			genes <- genes[1:min(top.n, length(genes))]
			cors.per.group <- lapply(split(sampleNames(object), pData(object)[[group]]), function(samples) {
				sapply(genes, function(gene) {
					cor(mirExprs[mir, samples], mrnaExprs[gene, samples], method='pearson', use='pairwise.complete.obs')
				})
			})

			boxplot(cors.per.group, main=paste(mir, '\n', paste(genes, collapse=', '), sep=''))
			return(cors.per.group)
		}
	})

	setMethod('correlationBoxplot', 'MatchedExpressionSet', function(object, mir, targets, group='molecular_subtype') {
		mirExprs  <- exprs(eset(object, 'mirna'))
		mrnaExprs <- exprs(eset(object, 'mrna'))

		mir <- match.arg(mir, rownames(mirExprs))

		if (is.null(targets)) {
			warning(mir, ' has no targets!')
		} else {

			cors.per.group <- lapply(split(sampleNames(object), pData(object)[[group]]), function(samples) {
				sapply(targets, function(gene) {
					cor(mirExprs[mir, samples], mrnaExprs[gene, samples], method='pearson', use='pairwise.complete.obs')
				})
			})

			boxplot(cors.per.group, main=paste(mir, '\n', paste(genes, collapse=', '), sep=''))
		}
	})
}

setMethod('plotProfile', 'MatchedExpressionSet', function(object, mir=featureNames(eset(object, 'mirna'))[1], sgID = sgIDs(object)[1], idfID = idfIDs(object)[1], top.n=10, group='molecular_subtype') {


	if (length(mir) > 1) {
		stop('You should provide only 1 mir!')
	}
#	mirExprs <- exprs(eset(object, 'mirna'))
#	mrnaExprs <- exprs(eset(object, 'mrna'))
	mirExprs <- exprs(filteredEset(object, 'mirna', sgID))
	mrnaExprs <- exprs(filteredEset(object, 'mrna', sgID))

	mir <- match.arg(mir, rownames(mirExprs))


	sample_ids <- sg(object, sgID)



	group_factor <- do.call('interaction', c(lapply(group, function(g) addNA(factor(eset(object, 'mirna')[[g]]))), sep=' - '))

	
	mirna_profile <- data.frame(
		sample_id  = sample_ids,
		group      = group_factor,
		expression = mirExprs[mir, ],
		mir        = mir,
		label      = mir,
		type       = 'mir'
	)



	interactions <- filteredIdf(object, idfID, sgID)
	interactions$cors <- correlations(object, sgID, idfID)
#	browser()
	interactions <- interactions[order(interactions$cors), ]

	genes <- interactions[[2]][interactions[[1]] == mir]
	if (length(genes) == 0) {
		stop('No interactions found for ', mir)
	}

 	genes <- genes[1:min(top.n, length(genes))]

 	genes.cors <- interactions$cors[match(genes, interactions[[2]])]

	mrna_profiles <- data.frame(
		sample_id  = rep(sample_ids, length(genes)),
		group      = rep(group_factor, length(genes)),
		expression = as.vector(t(mrnaExprs[genes, , drop=FALSE])),
		mir        = mir,
		label      = rep(paste(genes, ' (', signif(genes.cors, 5), ')', sep=''), each=length(sample_ids)),
		type       = paste('Top', top.n, 'targets')
	)

	profiles <- rbind(mirna_profile, mrna_profiles)

	## order the levels of 'sample_id' by group and expression
	profiles$sample_id <- factor(as.character(profiles$sample_id), levels=sample_ids[order(mirna_profile$group, mirna_profile$expression)])
	profiles$group <- factor(profiles$group)

#	profiles <- profiles[order(profiles$mir_rank),]



#	idx <- order(mirna_profile$expression)
#	profiles <- profiles[rep(idx, 1+length(genes)) + rep(0:length(genes), each=nrow(mirna_profile)) * nrow(mirna_profile), ]
#browser()


	# qplot(x=sample_id, y=expression, data=profiles, colour=factor(subtype))
#	g<-ggplot(aes(x=sample_id,y=expression,colour=subtype),data=profiles) + geom_point() + facet_wrap(~ mir + type, scale="free")#+geom_point(aes(),data=profiles_mRNA)
	g <- ggplot(aes(x=sample_id, y=expression, colour=label, group=label), data=profiles)
	g <- g + stat_smooth(method="lm", aes(fill = group, group=interaction(group, label)))
	g <- g + geom_point() + geom_line(aes(group=interaction(group, label)))
	g <- g + facet_wrap(~ type + mir, ncol=1, scale="free")#+geom_point(aes(),data=profiles_mRNA)
	g <- g + labs(title="miR interaction profile")
	g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))


	print(g)

	profiles
})

setMethod('plotEcdfs', 'MatchedExpressionSet', function(object, interaction.names=interactionNames(object), n.bins=1000, xlim = c(-.5, .5), alternative='less', make.pdf=FALSE, ...) {

	all.cors <- correlations(object) 

	breaks <- seq(-1, 1, length.out=n.bins+1)

	if (make.pdf) {
		d1 <- density(all.cors)
		ylim <- range(d1$y)
		ylim[2] <- ylim[2]*1.7
		plot(d1$x, d1$y,
			col  = 'black',
			type = 'l',
			lty = 2, lwd = 2,
			xlim = xlim,
			ylim = ylim,
			xlab = 'Correlation',
			ylab = 'Density',
			main = "Probability Density plot",
			...)

	} else {
		plot(breaks, ecdf(all.cors)(breaks),
			col  = 'black',
			type = 'l',
			lty = 2, lwd = 1,
			xlim = xlim,
			xlab = 'Correlation',
			ylab = 'Density',
			main = "Emperical Cumulative Density plot",
			...)
	}


	n <- length(interaction.names)
	names <- interaction.names

	shifts <- numeric(n)
	p.vals <- numeric(n)
	n.interactions <- numeric(n)

	for (i in 1:n) {
		cors.bg <- correlations(object, id=names[i], na.rm=TRUE, invert=TRUE) 
		cors.fg <- correlations(object, id=names[i], na.rm=TRUE, invert=FALSE)


		if (length(cors.fg) <= 1) {
			warning("Skipping '", names[i], "': to few correlations to the foreground make a density plot")
			next			
		}
		if (length(cors.bg) <= 1) {
			warning("Skipping '", names[i], "': to few correlations to the foreground make a density plot")
			next			
		}

		di <- density(cors.fg)

		if (make.pdf) {
			lines(di$x, di$y, lty=1, lwd=2, col=i+1)
		} else {
			lines(breaks, ecdf(cors.fg)(breaks), lty=1, lwd=1, col=i+1)
		}

		shifts[i] <- median(cors.fg) - median(cors.bg)
		p.vals[i] <- wilcox.test(cors.fg, cors.bg, alternative=alternative)$p.value
		n.interactions[i] <- length(cors.fg)
	}


	

	legend('topleft',
		paste(names, ' (n=', n.interactions,', shift=', signif(shifts, 3), ', p=', signif(p.vals, 3), ')', sep='')
		, lty=1, lwd=2, col=1+(1:n), cex=.8, bty='n'
	)

})




setMethod('plotCorrelations', 'MatchedExpressionSet', function(object, type=c('hist', 'ecdf', 'anders-fisher-z', 'anders-fraction'), n.bins=100, ...) {
	type <- match.arg(type, several.ok = TRUE)

	breaks <- seq(-1, 1, length.out=n.bins+1)

	for (name in interactionNames(object)) {

		cors.bg <- correlations(object, id=name, na.rm=TRUE, invert=TRUE) 
		cors.fg <- correlations(object, id=name, na.rm=TRUE, invert=FALSE) 


		n.bg <- length(cors.bg)
		n.fg <- length(cors.fg)

		if (length(cors.fg) == 0) {
			warning("Skipping '", name, "': all correlation values are NA.")
			next			
		}

		tt.less    <- t.test(cors.fg, cors.bg, alternative='less')
		tt.greater <- t.test(cors.fg, cors.bg, alternative='greater')

		wilcox.less    <- wilcox.test(cors.fg, cors.bg, alternative='less')
		wilcox.greater <- wilcox.test(cors.fg, cors.bg, alternative='greater')

		par(mfcol=n2mfrow(length(type) ))

		if ('hist' %in% type) {
			message('Plotting "hist" for: ', name)
			hist(cors.fg,
				border = 'darkblue',
				breaks = breaks,
				freq   = FALSE,
				xlab   = 'Correlation',
				main   = "Histogram"
			)
			lines(density(cors.bg), col='black', lty=2, lwd=1)
			legend('topright', c('background', name), lty=c(2, 1), lwd=c(1,3), col=c('black', 'darkblue'), bty='n')

			cm <- correlationMatrix(object)
			cm.n.str <- paste(paste(dim(cm), ' (', names(dimnames(cm)), ')' , sep='', collapse=' x '), ' = ', prod(dim(cm)), sep='')
			idf.n.str <- nrow(interactions(object, name))

			legend(
				#-1, 3,
				#adj=c(0, 1),
				'topleft',
				bty='n',
				cex=0.6,
				legend=c(
					'Data:',
					paste('- ', name(object)),
					paste('- ', esetNames(object)[1], ': ', annotation(eset(object, esetNames(object)[1])),  sep=''),
					paste('- ', esetNames(object)[2], ': ', annotation(eset(object, esetNames(object)[2])),  sep=''),
					'Interactions:',
					paste('- ID: ',   name,     sep=''),
					paste('- N: ',    idf.n.str, sep=''),
					'Correlation matrix:',
					paste('- N samples: ',     length(sampleNames(object)), sep=''),
					paste('- Dimensions: ',    cm.n.str, sep=''),
					'Correlations stats (after NA removal):',
					paste('- N foreground: ',        n.fg, sep=''),
					paste('- N background: ',        n.bg, sep=''),
					paste('- N total: ',             n.bg+n.fg, sep=''),
					paste('- Fraction: ',            signif(n.fg / (n.bg + n.fg), 3), sep=''),
					paste('- Ratio: ',               signif(n.fg / n.bg, 3), sep=''),
					'T-test:',
					paste('- Statistic: ',           signif(tt.less$statistic, 3), sep=''),
					paste('- P-value (less): ',      signif(tt.less$p.value, 3), sep=''),
					paste('- P-value (greater): ',   signif(tt.greater$p.value, 3), sep=''),
					'Wilcoxon-test:',
					paste('- Statistic: ',           signif(wilcox.less$statistic, 3), sep=''),
					paste('- P-value (less): ',      signif(wilcox.less$p.value, 3), sep=''),
					paste('- P-value (greater): ',   signif(wilcox.greater$p.value, 3), sep=''),
					'Distribution shift:',
					paste('- MED(bg) - MED(fg) ',      signif(median(cors.bg) - median(cors.fg), 3), sep=''),
					'Plot settings:',
					paste('- N bins: ',              n.bins, sep='')
				)
			)
		}


		if ('ecdf' %in% type) {
			message('Plotting "ecdf" for: ', name)

			plot(breaks, ecdf(cors.bg)(breaks),
				col  = 'black',
				type = 'l',
				lty = 2, lwd = 1,
				xlab = 'Correlation',
				ylab = 'Density',
				main = "Emperical Cumulative Density plot",
				...)
			
			lines(breaks, ecdf(cors.fg)(breaks), col='darkblue', lty=1, lwd=3)
		}


		if ('anders' %in% substring(type, 1, 6)) {
		
			bg.bins <- cut(cors.bg, breaks, include.lowest=TRUE)
			fg.bins <- cut(cors.fg, breaks, include.lowest=TRUE)

			bg.counts <- tapply(cors.bg, bg.bins, length)
			fg.counts <- tapply(cors.fg, fg.bins, length)

			fractions <- fg.counts/bg.counts
			na.idx <- is.na(fractions)

			fts <- sapply(1:n.bins, function(i) {
				if (na.idx[i]) {
					c(
						fisher.estimate = NA,
						fisher.pvalue   = NA
					)
				} else {
			
					ct <- matrix(c(
						fg.counts[i], bg.counts[i],
						n.fg - fg.counts[i],         n.bg - bg.counts[i]
					), nrow=2, dimnames = list(
                        target = c("yes", "no"),
						bin = c("in bin ", "not in bin")
                    ))

					ft <- fisher.test(ct, alternative='less')
					c(
						fisher.estimate = ft$estimate,
						fisher.pvalue   = ft$p.value
					)
				}
			})

			fisher.test.z.values <- qnorm(fts[2,])
			bin.means <- breaks[1:length(breaks)-1]+diff(breaks)/2

			if ('anders-fisher-z' %in% type) {
				message('Plotting "anders-fisher-z" for: ', name)
				plot(bin.means[!na.idx], fisher.test.z.values[!na.idx], type='l', lwd=3, col='darkblue', ylim=c(-10, 10), xlim=c(-1,1), main="Anders plot", xlab='Correlation', ylab='Fisher Test Z-score')
				abline(h=0, lty=2)
				abline(v=breaks, lty=3)
			}

			if ('anders-fraction' %in% type) {
				message('Plotting "anders-fraction" for: ', name)
				plot(bin.means[!na.idx], fractions[!na.idx], type='l', lwd=3, col='darkblue', xlim=c(-1,1), main="Anders plot", xlab='Correlation', ylab='Fraction')
				abline(h=length(cors.fg)/length(cors.bg), lty=2)
				abline(v=breaks, lty=3)
			}
			#legend('topright', c('background', name), lty=c(2, 1), lwd=c(1,3), col=c('black', 'darkblue'))
		}
	}
})







setMethod('wilcoxTest', 'MatchedExpressionSet', function(object, id, split=FALSE, alternative='less') {
	

	if (split) {

		idf <- interactions(object, id)
		splitlist <- split(idf[[2]], idf[[1]])

		c.mat <- correlationMatrix(object)

		# filter correlation matrix
		c.mat.melted <- melt(c.mat[unique(idf[[1]]), unique(idf[[2]]), drop=FALSE])

		wilcox.list <- lapply(names(splitlist), function(mir) {
			genes <- splitlist[[mir]]

			message('Wilcoxon test for: ', mir)

			idx <- c.mat.melted[[1]] %in% mir & c.mat.melted[[2]] %in% genes
			if (sum(idx) <= 1) {
				warning("0 or 1 correlation values for ", mir, " found!. Wilcoxon not possible.")
				NULL
			} else {
				cors.fg <- c.mat.melted$value[idx]
				cors.bg <- c.mat.melted$value[!idx]

				wilcox.test(cors.fg, cors.bg, alternative=alternative)
			}
		})

		names(wilcox.list) <- names(splitlist)
		wilcox.list

	} else {

		cors.bg <- correlations(object, id, na.rm=TRUE, invert=TRUE) 
		cors.fg <- correlations(object, id, na.rm=TRUE, invert=FALSE) 


		if (length(cors.fg) == 0) {
			warning("Skipping '", id, "': all correlation values are NA.")
			NULL
		} else {
			wilcox.test(x = cors.fg, y = cors.bg, alternative = alternative)
		}
	}
})







###
### Global Test methods
###

setMethod('globalTest', 'MatchedExpressionSet', function(object, genericSubsets=list()) {
	

	mir.exprs  <- exprs(eset(object, id = 'mirna'))
	gene.exprs <- exprs(eset(object, id = 'mrna'))

	mirs <- rownames(mir.exprs)
	genes <- rownames(gene.exprs)

	# filter genericSubsets
	genericSubsets <- lapply(genericSubsets, function(subset) {
		intersect(subset, genes)
	})
	genericSubsets[sapply(genericSubsets, function(x) length(x) == 0)] <- NULL


	specificSubsetsList <- lapply(object@interactionDataFrames, function(idf) {
		split(idf[[2]], idf[[1]])
	})

	
	gts <- lapply(mirs, function(mir) {
		
		specificSubsets <- lapply(specificSubsetsList, function(x) x[[mir]])
		
		subsets <- c(
			specificSubsets,
			genericSubsets
		)


		message('Global test for: ', mir, '...')

		# optimisation trick
		# filter gene.exprs for only the genes that are in the subset.

#		alt <- t(gene.exprs)
#		for (subset in subsets) {
#			if (length(subset) == 627) {
#				alt[, subset[100]] <- alt[, subset[100]] + rnorm(nrow(alt), 0, 1e-10)
#				message('Added some minor noise to avoid stalling of global test.')
#			}
#		}

#		lapply(names(subsets), function (subset_name) {
#			subset <- subsets[[subset_name]]
#			message(subset_name, ' n=', length(subset), '...')
#			gt(response=mir.exprs[mir, ], alternative=alt[, subset])			
#		})

#		res <- gt(response=mir.exprs[mir, ], alternative=alt, subsets=subsets)
		res <- gt(response=mir.exprs[mir, ], alternative=t(gene.exprs), subsets=subsets)
#		res <- gt(response=mir.exprs[mir, ], alternative=t(gene.exprs[unique(unlist(subsets)), ]), subsets=subsets)
		
		res@legend$cov <- c(
			paste("pos. assoc. with", mir),
			paste("neg. assoc. with", mir)
		)
		res@legend$subj <- c(
			paste("pos. residual", mir),
			paste("neg. residual", mir)
		)
		res
	})
	names(gts) <- mirs
	gts
})


setMethod('globalTestPheno', c('ExpressionSet'), function(object, interactionDataFrames, phenoColumns, nullPhenoColumns=list(), genericSubsets=list()) {
	#phenoColumn <- match.arg(phenoColumn, colnames(pData(object)))

	gene.exprs <- exprs(object)
	genes <- rownames(gene.exprs)

	X <- t(gene.exprs)

	genericSubsets <- lapply(genericSubsets, function(x) intersect(genes, x))

	genericSubsets <- c(genericSubsets, list(
		genes = genes
	))
	

	specificSubsetsList <- lapply(interactionDataFrames, function(idf) {
		split.df <- split(idf[[2]], idf[[1]])

		mapply(names(split.df), split.df, FUN=function(mir, genes) {
			genes
		})
	})

	specificSubsets <- do.call('c', specificSubsetsList)

	subsets <- c(
		genericSubsets,
		specificSubsets
	)


	results <- lapply(names(phenoColumns), function(phenoColumnName) {
		phenoColumn <- phenoColumns[[phenoColumnName]]
		nullPhenoColumns <- nullPhenoColumns[[phenoColumnName]]

		if (length(phenoColumn) == 1) {

			Y <- as.factor(pData(object)[[phenoColumn]])

			# remove NA's
			X <- X[!is.na(Y),]
			Y <- Y[!is.na(Y)]

			if (!is.null(nullPhenoColumns)) {
				null <- pData(object)[!is.na(Y), nullPhenoColumns]
			}

		} else if (length(phenoColumn) == 2) {
			# COX model
			# Assume right sensored survival data
			Y <- pData(object)[, phenoColumn]

			na.idx <- apply(Y, 1, function(x) any(is.na(x)))

			Y <- Y[!na.idx, ]
			X <- X[!na.idx, ]

			Y <- Surv(time=Y[[1]], event=Y[[2]])

			if (!is.null(nullPhenoColumns)) {
				null <- pData(object)[!na.idx, nullPhenoColumns]
			}

		}

		if ( length(Y) <= 2) {
			warning(paste(paste(phenoColumn, collapse='-'), ': 0 (non-NA)', sep=''))
			return(NULL)
		}

		
		if (!is.null(nullPhenoColumns)) {
			#null <- makeOrdinalFeatureMatrix(null)

			if (all(is.na(null))) {
				message(phenoColumnName)
				res <- gt(Y, X, subsets = subsets)
			} else {
				message(phenoColumnName, ', null: ', paste(nullPhenoColumns, collapse=', '))
				if (length(nullPhenoColumns) == 1) {
					if (is.character(null)) {
						#null <- sapply(unique(null), function(lev) as.integer(null == lev))
						null <- as.matrix(as.integer(factor(null)), ncol=1)
					}
				}
				res <- gt(Y, X, null=null, subsets = subsets)
			}
		} else {
			message(phenoColumnName)
			res <- gt(Y, X, subsets = subsets)
		}

		#res <- gt(Y, X, subsets = subsets)
		res@legend <- lapply(res@legend, function(x) gsub('Y', paste(phenoColumn, collapse='-'), x))
		res
	})
	names(results) <- names(phenoColumns)
	results
})



#
#
# TODO
# To fix the performance issue for LumA and the full data set
#
# filter the columns of both gene.exprs and mirna.exprs such that we only use mirna's and mrna genericSubsets
# Check also for individual mrna testing.
#
#
# Do phenoColumns, list of pheno columns
setMethod('globalTestPheno', 'MatchedExpressionSet', function(object, phenoColumns, nullPhenoColumns=list(), include.mirna=TRUE, genericSubsets=list()) {

	#phenoColumn <- match.arg(phenoColumn, colnames(pData(object)))


	gene.exprs <- exprs(eset(object, id = 'mrna'))
	genes <- rownames(gene.exprs)

	genericSubsets <- lapply(genericSubsets, function(x) intersect(genes, x))

	if (include.mirna) {
		mir.exprs  <- exprs(eset(object, id = 'mirna'))

		mirs <- rownames(mir.exprs)
		X <- t(rbind(mir.exprs, gene.exprs))

		genericSubsets <- c(genericSubsets, list(
			mirs_and_genes = c(mirs, genes),
			mirs = mirs,
			genes = genes
		))

	} else {
		X <- t(gene.exprs)

		genericSubsets <- c(genericSubsets, list(
			genes = genes
		))

	}

	specificSubsetsList <- lapply(object@interactionDataFrames, function(idf) {
		split.df <- split(idf[[2]], idf[[1]])

		mapply(names(split.df), split.df, FUN=function(mir, genes) {
			if(include.mirna) {
				c(mir, genes)
			} else {
				genes
			}
		})
	})

	specificSubsets <- do.call('c', specificSubsetsList)

	subsets <- c(
		genericSubsets,
		specificSubsets
	)


	results <- lapply(names(phenoColumns), function(phenoColumnName) {
		phenoColumn <- phenoColumns[[phenoColumnName]]

		nullPhenoColumns <- nullPhenoColumns[[phenoColumnName]]

		if (length(phenoColumn) == 1) {
			
			### FIX THIS!!
			#  Y <- as.factor(pData(object)[[phenoColumn]]) --> is w
			#
			####
			#Y <- as.factor(pData(object)[[phenoColumn]])
			Y <- pData(object)[[phenoColumn]]

			# remove NA's
			X <- X[!is.na(Y),]
			Y <- Y[!is.na(Y)]

			if (!is.null(nullPhenoColumns)) {
				null <- pData(object)[!is.na(Y), nullPhenoColumns]
			}


		} else if (length(phenoColumn) == 2) {
			# COX model
			# Assume right sensored survival data
			Y <- pData(object)[, phenoColumn]

			na.idx <- apply(Y, 1, function(x) any(is.na(x)))

			Y <- Y[!na.idx, ]
			X <- X[!na.idx, ]

			Y <- Surv(time=Y[[1]], event=Y[[2]])


			if (!is.null(nullPhenoColumns)) {
				null <- pData(object)[!na.idx, nullPhenoColumns]
			}
		}


		if ( length(Y) <= 2) {
			warning(paste(paste(phenoColumn, collapse='-'), ': 0 (non-NA)', sep=''))
			return(NULL)
		}

		if (!is.null(nullPhenoColumns)) {
			#null <- makeBinaryFeatureMatrix(null)
			null <- makeOrdinalFeatureMatrix(null)
			#browser()

			if (all(is.na(null))) {
				message(phenoColumnName)
				res <- gt(Y, X, subsets = subsets)
			} else {
				#browser()
				message(phenoColumnName, ', null: ', paste(nullPhenoColumns, collapse=', '))
				res <- gt(Y, X, null=null, subsets = subsets)
			}
		} else {
			message(phenoColumnName)
		#	if (phenoColumn == 'lymph_nodes_positive') {
		#		browser()
		#	}
			res <- gt(Y, X, subsets = subsets)
		}
		res@legend <- lapply(res@legend, function(x) gsub('Y', paste(phenoColumn, collapse='-'), x))
		res
	})
	names(results) <- names(phenoColumns)
	results
})


setMethod('thresholdVariationPlot', 'MatchedExpressionSet', function(
	object,
	id,
	fixed.eset               = c('mirna', 'mrna'),
	fixed.assayDataElement   = 'exprs',
	fixed.threshold          = NA,
	varying.assayDataElement = 'exprs',   # 'A'
	n.varying.thresholds     = 20,
  	varying.thresholds.log.scale = FALSE,
	alternative = 'less',
  	aggregation.function = mean,
	y = c('shift', 'p-value', 'n.varying'),
	...) {



	fixed.eset   <- match.arg(fixed.eset)
	y            <- match.arg(y, several.ok=TRUE)

	# opposite of the fixed.eset
	varying.eset <- switch(fixed.eset,
		mirna = 'mrna',
		mrna  = 'mirna'
	)

	# 
	fixed.means <- apply(assayDataElement(eset(mEset, fixed.eset), fixed.assayDataElement), 1, aggregation.function, na.rm=TRUE)

	if (!is.na(fixed.threshold)) {
		fixed.idx <- fixed.means >= fixed.threshold
		fixed.idx[is.na(fixed.idx)] <- FALSE
	} else {
		fixed.idx <- rep(TRUE, length(fixed.means))
	}

	
	varying.means <- apply(assayDataElement(eset(mEset, varying.eset), varying.assayDataElement), 1, aggregation.function, na.rm=TRUE)

	if (varying.thresholds.log.scale) {
		varying.thresholds <- 10^(seq(-7, -1, .25))

	} else {

		varying.thresholds <- seq(min(varying.means, na.rm=TRUE), max(varying.means, na.rm=TRUE), length.out=n.varying.thresholds)
		varying.thresholds <- varying.thresholds[1:(n.varying.thresholds-1)]
	}

	### Loop over thresholds

	result <- sapply(varying.thresholds, function(threshold) {
		varying.idx <- varying.means >= threshold
		varying.idx[is.na(varying.idx)] <- FALSE


		message('threshold ', threshold, ': ', sum(varying.idx), ' ', varying.eset, "'s...")


		subsetList <- list(fixed.idx, varying.idx)
		names(subsetList) <- c(fixed.eset, varying.eset)

		subsettedMatchedEset <- subsetByFeatures(mEset, subsetList)


		cors.bg <- correlations(subsettedMatchedEset, id, na.rm=TRUE, invert=TRUE) 
		cors.fg <- correlations(subsettedMatchedEset, id, na.rm=TRUE, invert=FALSE)
    
    
		if (length(cors.fg) < 2 | length(cors.fg) < 2) {
			c(
			shift          = NA,
			'p-value'      = NA,
			n.varying      = sum(varying.idx)
			)

		} else {

			p.val <- wilcox.test(cors.fg, cors.bg, alternative=alternative)$p.value
			c(
			shift          = median(cors.fg) - median(cors.bg),
			'p-value'      = -log10(p.val),
			n.varying      = sum(varying.idx)
			)
		}
	})

	main <- paste('Threshold varations - fixed ', fixed.eset, ' ', deparse(substitute(aggregation.function)), ' (', fixed.assayDataElement, ') at ', fixed.threshold, sep='')
	xlab <- paste('Threshold on ', varying.eset, ' ', deparse(substitute(aggregation.function)), ' (', varying.assayDataElement, ')', sep='')

	old.par <- par(mfrow=n2mfrow(length(y)))
	for (var in y) {
		if (varying.thresholds.log.scale) {
			plot(varying.thresholds, result[var, ], log='x', xaxt='n', type='l', main=main, xlab=xlab, ylab=var)
			axis(1, c(.1, .01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001))
		} else {
			plot(varying.thresholds, result[var, ], type='l', main=main, xlab=xlab, ylab=var)
		}
	}
	par(old.par)

})


