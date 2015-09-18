###
### Helper methods
###


## Rank test from the Endeavour paper (http://www.nature.com/nbt/journal/v24/n5/full/nbt1203.html)
## r: Vector of rank ratios. A rank ratio is the rank divided by the total number of ranked items (value between 0 and 1)
qrank <- function(r, verbose=FALSE) {

	r <- r[!is.na(r)]      # remove NA's

	if (length(r) == 0) {
		1
	} else if (length(r) == 1) {
		r
	} else {

		r <- sort(r)
		d <- length(r)

		if (verbose) {
			cat('----------------\n')
			cat("r =", r, '\n')
			cat("d =", d, '\n')
			cat('----------------\n\n')
		}

		w <- numeric(length=d + 1)
		w[1] <- 1
		w[2] <- r[d]

		for (k in 3:(d+1)) {
			w[k] <- 0
			f <- -1
			for (j in 0:(k-2)) {
				f = -(f * (k - j - 1) * r[d - k + 2]) / (j + 1)
				w[k] = w[k] + w[k - j - 1] * f

				if (verbose) {
					cat("k =", k, '\n')
					cat("j =", j, '\n')
					cat("f =", f, '\n')
					cat("w =", w, '\n\n')
				}
			}
		}

		w[d + 1]
	}
}



## Rank test for a data frame of scores
rank.test <- function(score.df, return.rank.ratio.matrix=FALSE, n.permutation=10000, verbose=TRUE) {
	
	# Compute the Rank ration matrix
	rank.ratio.matrix <- sapply(score.df, function(x) rank(x) / nrow(score.df))

	# and get the qvalues
	qvals <- apply(rank.ratio.matrix, 1, qrank)

	# do permutation to get pvalues
	if (verbose) {
		cat(n.permutation, ' permutations...\n')
	}
	qvals.null <- sapply(1:n.permutation, function(i) qrank(apply(rank.ratio.matrix, 2, sample, size=1)))
	pvals      <- sapply(qvals, function(q) sum(qvals.null < q) / length(qvals.null))
	pvals.bh   <- p.adjust(pvals, method='fdr')

	na.idx <- apply(as.matrix(score.df), 1, function(x) any(is.na(x)))
	rank.ratio.matrix[na.idx, ] <- NA
	qvals[na.idx] <- NA
	pvals[na.idx] <- NA
	pvals.bh[na.idx] <- NA

	if (return.rank.ratio.matrix) {
		mat <- cbind(
			rank.ratio.matrix,
			qvals,
			pvals,
			pvals.bh
		)
		colnames(mat) <- c(
			paste('rank-ratio-', colnames(score.df), sep=''),
			'rank-test-q-value',
			'rank-test-p-value',
			'rank-test-p-value-fdr'
		)
		mat
	} else {
		pvals
	}
}





## Load the parclip targets
# Deprecated
#loadPC <- function() {
#	parclip_targets <- read.csv(file='data/parclip/PAR_CLIP.csv', stringsAsFactors=FALSE)
#	parclip_targets <- split(parclip_targets$GeneID, parclip_targets$miRNA)
#	parclip_targets <- lapply(parclip_targets, unique)
#
#	mapping_file <- 'data/targetscan/miR_Family_Info.txt'
#	mapping <- read.delim(file=mapping_file, stringsAsFactors=FALSE)
#	fam2mirbase <- split(mapping$MiRBase.ID, mapping$miR.family)
#
#	parclip_targets <- lapply(fam2mirbase, function(x) {
#		y <- unique(unlist(parclip_targets[x]))
#		y[y != "" & y != "-"]
#	})
#	parclip_targets[sapply(parclip_targets, is.null) ] <- NULL
#	parclip_targets
#}


## Get the interaction for a model
getModelInteractions <- function(matched_annotation, fIDs, sgID, idfID) {

	regression_method <- strsplit(idfID, '-')[[1]][[2]]
	featureSet        <- strsplit(idfID, '-')[[1]][[3]]
	cutoff            <- as.numeric(strsplit(idfID, '-')[[1]][[4]])
	
	model_dir       <- outputDir(relative_output_dir=modelDir(regression_method, featureSet), matched_annotation, fIDs, sgID)
	prediction_file <- paste(model_dir, 'prediction.csv', sep='')
	
	preds <- read.csv(file=prediction_file, stringsAsFactors=FALSE)
	preds <- preds[!is.na(preds$prediction), ]
	
	# apply threshold
	preds <- unique(preds[preds$prediction > cutoff, c(1,2)])
	preds
}



## Add model interacion dataframes to mEset for model prediction
addModelPredictionIdfForValidation <- function(mEset, regression_method, featureSet, annotation, sgID, threshold, fIDs, individuals=TRUE) {
	
	if (individuals) {

		# split targetscan-all
		split.ts <- split(mEset@interactionDataFrames[['targetscan_predicted_conserved-all']], mEset@interactionDataFrames[['targetscan_predicted_conserved-all']]$miR.Family)

		for (mirna in names(split.ts)) {
			mEset@interactionDataFrames[[paste('targetscan_predicted_conserved', mirna, sep='_')]] <- NULL
			setInteractions(mEset, idfID=paste('targetscan_predicted_conserved', mirna, sep='_'), interactions=split.ts[[mirna]])
		}


		# split parclip-all
		split.pc <- split(mEset@interactionDataFrames[['parclip-all']], mEset@interactionDataFrames[['parclip-all']]$miR.Family)

		for (mirna in names(split.pc)) {
			mEset@interactionDataFrames[[paste('parclip', mirna, sep='_')]] <- NULL
			setInteractions(mEset, idfID=paste('parclip', mirna, sep='_'), interactions=split.pc[[mirna]])
		}
	}


	model_dir       <- paste("results/models/", regression_method,"_model-",featureSet, "/", annotation, '/', sgID, paste(rep('.', length(fIDs)), fIDs, collapse='', sep=''), '/', sep='')

	prediction_file <- paste(model_dir, 'prediction.csv', sep='')
	

	idfID <- paste('model_prediction', threshold, sep='_')
	mEset@interactionDataFrames[[idfID]] <- NULL


	preds <- read.csv(file=prediction_file, stringsAsFactors=FALSE)
	preds <- preds[!is.na(preds$prediction), ]
	preds <- preds[preds[[1]] %in% featureNames(mEset@expressionSets[[1]]), ]
	preds <- preds[preds[[2]] %in% featureNames(mEset@expressionSets[[2]]), ]
	# apply threshold
	preds <- unique(preds[preds$prediction > threshold, c(1,2)])


	setInteractions(mEset, idfID=idfID, interactions=preds)

	if (individuals) {
		# split model predictions
		split.preds <- split(preds, preds$miRNA.family)

		for (mirna in names(split.preds)) {
			mEset@interactionDataFrames[[paste(idfID, mirna, sep='_')]] <- NULL
			setInteractions(mEset, idfID=paste(idfID, mirna, sep='_'), interactions=split.preds[[mirna]])
		}
	}

	mEset
}


plotPredictionMeasures <- function(predictions, labels, cutoff) {
	## prediction on the 
	prediction_object <- prediction(predictions, labels)


	## performances and ideal performance measures
	tps <- sapply(predictions, function(t) sum( (predictions >= t) & labels))
	fps <- sapply(predictions, function(t) sum( (predictions >= t) & !labels))
	predictions_order_idx <- order(predictions)

	ideal_tps <- c(rep(sum(labels), sum(!labels)), sum(labels):1)
	ideal_fps <- c(sum(!labels):1, rep(0, sum(labels)))


	## plots

	# TP
	plot(predictions[predictions_order_idx], tps[predictions_order_idx], xlab='Cutoff', ylab='# True positives', type='l')
	lines(predictions[predictions_order_idx], ideal_tps, col='red')
	abline(v=cutoff, lty=2)

	# FP
	plot(predictions[predictions_order_idx], fps[predictions_order_idx], xlab='Cutoff', ylab='# False positives', type='l')
	lines(predictions[predictions_order_idx], ideal_fps, col='red')
	abline(v=cutoff, lty=2)

	# ROC
	plot(performance(prediction_object, "tpr", "fpr"))
	legend('bottomright', paste('AUC =', signif(performance(prediction_object, "auc")@y.values[[1]], 4)))

	# TPR
	plot(performance(prediction_object , "tpr"))
	abline(v=cutoff, lty=2)

	# FPR
	plot(performance(prediction_object, "fpr"))
	abline(v=cutoff, lty=2)

	# PPV
	plot(performance(prediction_object, "ppv"))
	abline(v=cutoff, lty=2)

	prediction_object
}





makeBinaryFeatureMatrix <- function(df) {
	bin.matrix <- do.call('cbind', lapply(df, function(col) {
		levels <- levels(factor(col))
		if (all(is.na(col))) {
			NULL
		} else if (length(levels) == 1) {
			NULL
		} else {
			#if (length(levels) == 2) {
			#	numeric.mat <- matrix(as.numeric(factor(col)) - 1, ncol=1)
			#	colnames(numeric.mat) <- levels[2]
			#} else {
				numeric.mat <- do.call('cbind', lapply(levels, function(level) {
					as.numeric(col == level)
				}))
				colnames(numeric.mat) <- levels
				#numeric.mat[,]
			#}
			numeric.mat[,-(ncol(numeric.mat)-1)]
		}
	}))
	rownames(bin.matrix) <- rownames(df)
	bin.matrix
}




makeModelMatrix <- function (formula, data, subset, weights, na.action, method = "qr", 
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...) 
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame") 
        return(mf)
    else if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
            3) else numeric(), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
            0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }
    }
    else {
    	browser()
        x <- model.matrix(mt, mf, contrasts)
        return (x)

        z <- if (is.null(w)) 
            lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
                ...)
        else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
            ...)
    }
    class(z) <- c(if (is.matrix(y)) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- x
    if (ret.y) 
        z$y <- y
    if (!qr) 
        z$qr <- NULL
    z
}


