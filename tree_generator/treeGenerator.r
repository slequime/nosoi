treeGenerator = function(transmission_matrix, samples, sampling_times)
	{
		output_format="tree"; nberOfTraits = dim(transmission_matrix)[2]-3	; trait_names = c()
		for (i in 4:dim(transmission_matrix)[2]) trait_names = c(trait_names, colnames(transmission_matrix)[i])
		nodes = samples; times = sampling_times; clades = nodes
		for (i in 1:length(nodes))
			{
				clade = paste0(nodes[i],"[&ind=",nodes[i])
				for (j in 1:nberOfTraits)
					{
						index = which(transmission_matrix[,"ind.ID"]==nodes[i])
						clade = paste0(clade,",",trait_names[j],"=",transmission_matrix[index,trait_names[j]])
					}
				clades[i] = paste0(clade,"]")
			}
		t0 = max(sampling_times)
		for (t in t0:0)
			{
				newNodes = c(); newTimes = c(); newClades = c()
				indices1 = rep(NA, length(nodes)); infDates = rep(NA, length(nodes)); ancestors = rep(NA, length(nodes))
				for (i in 1:length(nodes))
					{
						indices1[i] = which(transmission_matrix[,"ind.ID"]==nodes[i])
						infDates[i] = transmission_matrix[indices1[i], "inf.date"]
						ancestors[i] = transmission_matrix[indices1[i], "inf.by"]
					}
				indices2 = which(infDates==t)
				for (i in 1:length(indices2))
					{
						nodes[indices2[i]] = ancestors[indices2[i]]
					}
				indices_already_used = c()
				for (i in 1:length(nodes))
					{
						if (!i%in%indices_already_used)
							{
								indices3 = which(nodes==nodes[i])
								indices3 = indices3[which(indices3!=i)]
								if (length(indices3) > 0)
									{
										t1 = times[i]-t
										clade = paste0("(",clades[[i]],":",t1)
										for (j in c(indices3))
											{
												t2 = times[j]-t
												clade = paste0(clade,",",clades[[j]],":",t2)
											}
										clade = paste0(clade,")")
										ancestor = nodes[i]
										traits = paste0("[&ind=",ancestor)
										for (j in 1:nberOfTraits)
											{
												index = which(transmission_matrix[,"ind.ID"]==nodes[i])
												traits = paste0(traits,",",trait_names[j],"=",transmission_matrix[index,trait_names[j]])
											}
										clade = paste0(clade,traits,"]")
										newNodes = c(newNodes, nodes[i])
										newTimes = c(newTimes, t)
										newClades = c(newClades, clade)
										indices_already_used = c(indices_already_used, i, indices3)
									}	else		{
										newNodes = c(newNodes, nodes[i])
										newTimes = c(newTimes, times[i])
										newClades = c(newClades, clades[i])
									}
							}
					}
				nodes = newNodes; times = newTimes; clades = newClades; # print(length(nodes))
			}
		if (output_format == "tree") tree = read.tree(text=paste0(clades[[1]],";"))
		if (output_format == "text") tree = paste0(clades[[1]],";")
		return(tree)
	}
