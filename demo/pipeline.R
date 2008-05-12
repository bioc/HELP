	require("Biobase")
	require("HELP")
	if (toupper(readline("Run HELP pipeline (y/n)? ")) != "Y") {
		stop()
	}
	if (!is.null(readline(" ... hit <Return> to continue ... "))) {
		cat("Reading data files ... ")
		data.path <- system.file("pipeline.data", package="HELP")
		msp1.files <- dir(path=data.path, pattern="532[._]pair", full.names=TRUE)
		hpa2.files <- dir(path=data.path, pattern="635[._]pair", full.names=TRUE)
		N <- length(msp1.files)
		for (i in 1:N) {
			if (i == 1) {
				pairs <- readPairs(msp1.files[i],hpa2.files[i])
			}
			else {
				pairs <- readPairs(msp1.files[i],hpa2.files[i],pairs)
			}
		}
		chips <- sub("[_.]*532[_.]*pair.*", "", basename(msp1.files))
		samples.file <- file.path(path=data.path, file="sample.key.txt")
		chips <- readSampleKey(file=samples.file, chips=chips)
		sampleNames(pairs) <- chips
		ngd.file <- file.path(path=data.path, file="HELP.ngd.txt")
		ndf.file <- file.path(path=data.path, file="HELP.ndf.txt")
		pairs <- readDesign(ndf.file, ngd.file, pairs)
		exprs(pairs) <- log2(exprs(pairs))
		exprs2(pairs) <- log2(exprs2(pairs))
		cat("FINISHED\n")
	}
	if (toupper(readline("Calculate Tm values for probe set (y/n)? ")) == "Y") {
		cat("Calculating Tm values for ", dim(pairs)[1], " probes ... ", sep="")
		pairs <- calcTm(pairs)
		plot(density(as.numeric(getFeatures(pairs, "TM"))), main="Melting temperature distribution", xlab="Tm (Celcius)")
		dev.print(pdf,file="DEMO.Tm.density.pdf")
		cat("FINISHED\n")
	}
	if (toupper(readline("Calculate %GC values for probe set (y/n)? ")) == "Y") {
		cat("Calculating %GC values for ", dim(pairs)[1], " probes ... ", sep="")
		pairs <- calcGC(pairs)
		plot(density(as.numeric(getFeatures(pairs, "GC"))), main="GC-content distribution", xlab="GC-content (%)")
		dev.print(pdf,file="DEMO.GC.density.pdf")
		cat("FINISHED\n")
	}
	if (!is.null(readline(" ... hit <Return> to continue ... "))) {
		cat("Calculating cutoffs ... ")
		meth.cutoffs <- failed.cutoffs <- failed2.cutoffs <- rep(0,N)
		rand <- which(getFeatures(pairs, "TYPE") == "RAND")
		for (i in 1:N) {
			hpa <- getSamples(pairs, i, element="exprs2")
			msp <- getSamples(pairs, i, element="exprs")
			msp.rand <- msp[rand]
			hpa.rand <- hpa[rand]
			failed.cutoffs[i] <- median(msp.rand) + 2.5*mad(msp.rand)
			failed2.cutoffs[i] <- median(hpa.rand) + 2.5*mad(hpa.rand)
			meth.cutoffs[i] <- max(quantile(hpa.rand, 0.99), failed2.cutoffs[i])
		}
		cat("FINISHED\n")
		temp <- rbind(failed.cutoffs, failed2.cutoffs, meth.cutoffs)
		colnames(temp) <- chips
		rownames(temp) <- c("MspI (failed)", "HpaII (failed)", "HpaII (meth)")
		print(temp)
		rm(temp)
	}
	if (toupper(readline("Perform quality assessment of raw data (y/n)? ")) == "Y") {
		cat("Calculating prototypes ... ")
		msp.prototype <- calcPrototype(pairs, element="exprs")
		hpa.prototype <- calcPrototype(pairs, element="exprs2")
		ratio.prototype <- calcPrototype(getSamples(pairs, element="exprs2") - getSamples(pairs, element="exprs"))
		cat("FINISHED\n")
		cat("Plotting chip reconstruction for MspI prototype ... ")
		plotChip(pairs, msp.prototype, range=range(msp.prototype))
		dev.print(pdf, file="DEMO.prototype.image.pdf")
		cat("FINISHED\n")
		cat("Plotting chip reconstruction for sample 'Brain2' (versus MspI prototype) ... ")
		msp.Brain2 <- getSamples(pairs, "Brain2", element="exprs")
		plotChip(pairs, msp.Brain2 - msp.prototype)
		dev.print(pdf, file="DEMO.chip.image.pdf")
		cat("FINISHED\n")
	}
	if (toupper(readline("Plot fragment size v. signal intensity (y/n)? ")) == "Y") {
		cat("Generating plot for 'Brain2' ... ")
		plotFeature(pairs, sample="Brain2")
		dev.print(pdf,file="DEMO.size.v.intensity.pdf")
		cat("FINISHED\n")
	}
	pairs.norm <- pairs
	if (toupper(readline("Normalize raw data for chip set (y/n)? ")) == "Y") {
		cat("Normalizing data ... ")
		norand <- (getFeatures(pairs, "TYPE") == "DATA")
		size <- as.numeric(getFeatures(pairs, "SIZE"))
		for (i in 1:N) {
			hpa <- getSamples(pairs, i, element="exprs2")
			msp <- getSamples(pairs, i, element="exprs")
			nofail <- which(norand & (msp>failed.cutoffs[i] | hpa>failed2.cutoffs[i]))
			exprs(pairs.norm)[nofail, i] <- quantileNormalize(msp[nofail], size[nofail])
			meth <- which(norand & msp>failed.cutoffs[i] & hpa<=meth.cutoffs[i])
			exprs2(pairs.norm)[meth, i] <- quantileNormalize(hpa[meth], size[meth])
			nometh <- which(norand & hpa>meth.cutoffs[i])
			exprs2(pairs.norm)[nometh, i] <- quantileNormalize(hpa[nometh], size[nometh])
		}
		cat("FINISHED\n")
	}
	if (toupper(readline("Show windowing scheme and data distribution per bin (y/n)? ")) == "Y") {
		cat("Plotting bins ... ")
		size <- as.numeric(getFeatures(pairs, "SIZE"))
		plotBins(pairs,size,sample="Brain1",element="exprs2")
		dev.print(pdf,file="DEMO.bins.pdf")
		cat("FINISHED\n")
	}
	seqids <- getFeatures(pairs.norm, "PROBE_ID")
	ratios.combined <- getSamples(pairs.norm, element="exprs2") - getSamples(pairs.norm, element="exprs")
	rownames(ratios.combined) <- seqids
	if (toupper(readline("Combine probe-level data by HpaII fragment (y/n)? ")) == "Y") {
		cat("Combining data for chip set ... ")
		msp.combined <- combineData(pairs.norm, seqids, element="exprs")
		hpa.combined <- combineData(pairs.norm, seqids, element="exprs2")
		ratios.combined <- combineData(getSamples(pairs.norm, element="exprs2") - getSamples(pairs.norm, element="exprs"), seqids)
		cat("FINISHED\n")
	}
	if (toupper(readline("Plot pairwise comparison of arrays in dataset (y/n)? ")) == "Y") {
		cat("Plotting array pairs ... ")
		plotPairs(ratios.combined)
		cat("FINISHED\n")
		dev.print(pdf,file="DEMO.tree.pairs.pdf")
	}
	if (toupper(readline("Create wiggle tracks for data (y/n)? ")) == "Y") {
		cat("Writing wiggle tracks  ... ")
		map <- match(rownames(ratios.combined), seqids)
		chr <- getFeatures(pairs, "CHR")[map]
		start <- getFeatures(pairs, "START")[map]
		end <- getFeatures(pairs, "STOP")[map]
		createWiggle(ratios.combined, cbind(chr, start, end), file="./DEMO.ratios.wig")
		cat("FINISHED\n")
		rm(map, chr, start, end)
	}
