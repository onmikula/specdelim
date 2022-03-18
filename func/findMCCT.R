findMCCT <- function(trees, file=NULL, return=TRUE, digits=2, method="mean", monophyletic=FALSE) {

	print("findMCCT is deprecated, use find_mcct at https://github.com/onmikula/find_mcct instead")

	mcct <- phangorn::maxCladeCred(trees, tree=TRUE, part=NULL, rooted=TRUE)
	if (method != "topology" & !is.null(method)) {
		mcct <- comm_anc_heights(phy=mcct, trees=trees, method=method, monophyletic=monophyletic)
	}
	mcct$node.label <- formatC(probClades(phy=mcct, psample=trees), format="f", digits=digits)
	if (!is.null(file)) {
		ape::write.tree(mcct, file)
	}
	if (isTRUE(return)) {
		return(mcct)
	}
}


comm_anc_heights <- function(phy, trees, monophyletic=FALSE, method="mean") {
	fun <- match.fun(method)
	R <- Ntip(phy) + 1
	tips <- phy$tip.label
	torders <- lapply(lapply(trees, "[[", "tip.label"), order)
	theights <- matrix(, length(tips), length(trees))
	for (i in seq_along(trees)) {
		tpaths <- ape::nodepath(trees[[i]])
		tbranches <- lapply(tpaths, function(x) match(x, trees[[i]]$edge[,2]))
		theights[,i] <- sapply(tbranches, function(x) sum(trees[[i]]$edge.length[x], na.rm=TRUE))[torders[[i]]]
	}
	TH <- apply(theights, 1, fun)
	TH <- TH[match(tips, sort(tips))]
	clades <- lapply(prop.part(phy), function(x) tips[x])
	nheights <- matrix(, length(clades), length(trees))
	for (i in seq_along(trees)) {
		if (isTRUE(monophyletic)) {
			monophyl <- which(sapply(clades, function(x) ape::is.monophyletic(trees[[i]], x)))
		} else {
			monophyl <- rep_len(TRUE, length(clades))
		}
		anc <- sapply(clades[monophyl], function(x) getMRCA(trees[[i]], x))
		npaths <- lapply(anc, function(a) ape::nodepath(trees[[i]], from=R, to=a))
		nbranches <- lapply(npaths, function(x) match(x, trees[[i]]$edge[,2]))
		nheights[monophyl,i] <- sapply(nbranches, function(x) sum(trees[[i]]$edge.length[x], na.rm=TRUE))
	}
	NH <- apply(nheights, 1, fun, na.rm=TRUE)
	H <- c(TH, NH)
	phy$edge.length <- as.numeric(diff(t(matrix(H[phy$edge], Nedge(phy), 2))))
	return(phy)
}


probClades <- function(phy, psample, complete=FALSE, species=NULL) {
	clades <- ape::prop.part(phy)
	clades <- lapply(clades, function(x) attr(clades, "labels")[x])
	clades <- lapply(clades, sort)
	part <- ape::prop.part(psample)
	pppart <- attr(part, "number") / length(psample)
	part <- lapply(part, function(x) attr(part, "labels")[x])
	part <- lapply(part, sort)
	ord <- sapply(seq_along(clades),  function(i) which(sapply(part, identical, y=clades[[i]])))
	pp <- pppart[ord]
	if (isTRUE(complete)) {
		if (!is.null(species)) {
			binarize <- function(x, maxn) as.numeric(seq(maxn) %in% x)
			part <- lapply(lapply(part, function(x) phy$tip.label[x]), match, table=species)
			part <- apply(sapply(part, binarize, maxn=length(species)), 2, paste, collapse="")
		}
		attr(pp, "part") <- data.frame(part=part, pp=ppart, stringsAsFactors=FALSE)
	}
	return(pp)
}



