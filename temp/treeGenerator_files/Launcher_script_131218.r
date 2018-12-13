library(ape)
library(phytools)
library(RColorBrewer)
library(OutbreakTools)

source("treeGenerator.r")

transmission_matrix = read.csv("transmission.csv", header=T)
transmission_matrix$hosts.ID = as.character(transmission_matrix$hosts.ID)
transmission_matrix$inf.by = as.character(transmission_matrix$inf.by)

nberOfRandomSamples = 100
samples_indices = sample(1:dim(transmission_matrix)[1], nberOfRandomSamples)
samples = as.character(transmission_matrix[samples_indices,"hosts.ID"])
last_infection_times = rep(NA, length(samples))
for (i in 1:length(samples))
	{
		indices = which(transmission_matrix[,"inf.by"]==samples[i])
		if (length(indices) > 0)
			{
				last_infection_times[i] = max(transmission_matrix[indices,"inf.time"])
			}	else		{
				last_infection_times[i] = transmission_matrix[samples_indices[i],"inf.time"]
			}
	}
sampling_times = last_infection_times + 1

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
axis(1, pos=-1, lwd.tick=0.5, cex.axis=0.6, lwd=0.5, tck=-0.005, col.axis="gray30", mgp=c(0,-0.1,0)); dev.off()
# axisPhylo(pos=-1, lwd.tick=0.5, cex.axis=0.6, lwd=0.5, tck=-0.005, col.axis="gray30", mgp=c(0,-0.1,0)); dev.off()
