library(ape)
library(phytools)
library(RColorBrewer)
library(OutbreakTools)

source("treeGenerator.r")

humans = read.csv("humans.csv", header=T); vectors = read.csv("vectors.csv", header=T)

transmission_matrix = data.frame(matrix(nrow=dim(humans)[1]+dim(vectors)[1], ncol=4))
	# the three first columns of the matrix are compulsory, the next ones are annotations
colnames(transmission_matrix) = c("ind.ID","inf.by","inf.date","inf.loc")
transmission_matrix$ind.ID = c(as.character(humans$hosts.ID), as.character(vectors$ID.vect))
transmission_matrix$inf.by = c(as.character(humans$inf.by), as.character(vectors$inf.by))
transmission_matrix$inf.date = c(as.numeric(humans$inf.date), as.numeric(vectors$infected.date))
transmission_matrix$inf.loc = c(as.character(humans$inf.in), as.character(vectors$loc))

nberOfRandomSamples = 100
samples_indices = sample(1:dim(transmission_matrix)[1], nberOfRandomSamples)
samples = transmission_matrix[samples_indices,"ind.ID"]
sampling_times = transmission_matrix[samples_indices,"inf.date"]
for (i in 1:length(sampling_times))
	{
		sampling_times[i] = sampling_times[i]+sample(c(5:10),1)
	}

tree = treeGenerator(transmission_matrix, samples, sampling_times)
write.tree(tree, "phylogeny.tree")

pdf("phylogeny.pdf"); col0 = colorRampPalette(brewer.pal(9,'PuBu'))(101)[1]
cols = colorRampPalette(brewer.pal(9,'PuBu'))(101)[(((nodeHeights(tree)[,2])/(max(nodeHeights(tree))-min(nodeHeights(tree))))*100)+1]
par(oma=c(1.3,0.5,0,0.5), mar=c(0,0,0,0), mgp=c(1,0.1,0)); plot(tree, show.tip.label=T, edge.width=0.5, cex=0.4, align.tip.label=3)
for (j in 1:dim(tree$edge)[1])
	{
		if (j == 1)
			{
				nodelabels(node=tree$edge[j,1], pch=16, cex=0.9, col=col0)
				nodelabels(node=tree$edge[j,1], pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
			}
		nodelabels(node=tree$edge[j,2], pch=16, cex=0.9, col=cols[j])
		nodelabels(node=tree$edge[j,2], pch=1, cex=0.9, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
	}
axisPhylo(pos=-1, lwd.tick=0.5, cex.axis=0.6, lwd=0.5, tck=-0.005, col.axis="gray30", mgp=c(0,-0.1,0)); dev.off()
