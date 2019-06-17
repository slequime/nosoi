# treeGenerator = function(transmission_matrix, samples, sampling_times)
# {
#   output_format = "tree"; nberOfTraits = dim(transmission_matrix)[2]-3; trait_names = c()
#   for (i in 4:dim(transmission_matrix)[2]) trait_names = c(trait_names, colnames(transmission_matrix)[i])
#   nodes = samples; times = sampling_times; clades = nodes
#   for (i in 1:length(nodes))
#   {
#     clade = paste0(nodes[i],"[&ind=",nodes[i])
#     for (j in 1:nberOfTraits)
#     {
#       index = which(transmission_matrix[,"hosts.ID"]==nodes[i])
#       clade = paste0(clade,",",trait_names[j],"=",transmission_matrix[index,trait_names[j]])
#     }
#     clades[i] = paste0(as.character(clade),"]")
#   }
#   t0 = max(sampling_times)
#   for (t in t0:0)
#   {
#     newNodes = c(); newTimes = c(); newClades = c()
#     indices1 = rep(NA, length(nodes)); infDates = rep(NA, length(nodes)); ancestors = rep(NA, length(nodes))
#     for (i in 1:length(nodes))
#     {
#       indices1[i] = which(transmission_matrix[,"hosts.ID"]==nodes[i])
#       infDates[i] = transmission_matrix[indices1[i], "inf.time"]
#       ancestors[i] = transmission_matrix[indices1[i], "inf.by"]
#     }
#     indices2 = which(infDates==t)
#     for (i in 1:length(indices2))
#     {
#       nodes[indices2[i]] = ancestors[indices2[i]]
#     }
#     indices_already_used = c()
#     for (i in 1:length(nodes))
#     {
#       if (!i%in%indices_already_used)
#       {
#         indices3 = which(nodes==nodes[i])
#         indices3 = indices3[which(indices3!=i)]
#         if (length(indices3) > 0)
#         {
#           t1 = times[i]-t
#           clade = paste0("(",clades[[i]],":",t1)
#           for (j in c(indices3))
#           {
#             t2 = times[j]-t
#             clade = paste0(clade,",",clades[[j]],":",t2)
#           }
#           clade = paste0(clade,")")
#           ancestor = nodes[i]
#           traits = paste0("[&ind=",ancestor)
#           for (j in 1:nberOfTraits)
#           {
#             index = which(transmission_matrix[,"hosts.ID"]==nodes[i])
#             traits = paste0(traits,",",trait_names[j],"=",transmission_matrix[index,trait_names[j]])
#           }
#           clade = paste0(clade,traits,"]")
#           newNodes = c(newNodes, nodes[i])
#           newTimes = c(newTimes, t)
#           newClades = c(newClades, clade)
#           indices_already_used = c(indices_already_used, i, indices3)
#         }	else		{
#           newNodes = c(newNodes, nodes[i])
#           newTimes = c(newTimes, times[i])
#           newClades = c(newClades, clades[i])
#         }
#       }
#     }
#     nodes = newNodes; times = newTimes; clades = newClades; print(length(nodes))
#   }
#   if (output_format == "tree") tree = read.tree(text=paste0(clades[[1]],";"))
#   if (output_format == "text") tree = paste0(clades[[1]],";")
#   return(tree)
# }

getTransmissionTree <- function(nosoiInf) {
  if (!requireNamespace("ape", quietly = TRUE) || !requireNamespace("tidytree", quietly = TRUE) || !requireNamespace("treeio", quietly = TRUE)) {
    stop("Packages 'ape', 'tidytree' and 'treeio' are needed for transmission tree generation.",
         call. = FALSE)
  }

  table.hosts <- nosoiInf$host.info.A$table.hosts
  setorder(table.hosts, "inf.time")

  # Indicators for tips and nodes
  # each host gives a tip
  table.hosts[, "indTips" := .I]
  # each transmission event gives a node
  table.hosts[, "indNodes" := .GRP, by = c("inf.by", "inf.time")]
  # if NA then last till the end
  table.hosts[is.na(get("out.time")), ("out.time") := nosoiInf$total.time]
  # Names of the hosts
  hosts <- table.hosts[["hosts.ID"]]

  # Characteristics of the tree
  nTips <- nrow(table.hosts)
  nNode <- length(unique(table.hosts[["label"]])) - 1
  nEdges <- nNode + nTips - 1

  # Initialize the tree
  treeTable <- tidytree::tibble(parent = NA_integer_, node = NA_integer_, branch.length = NA,
                                label = NA_character_, host = NA, state = NA,
                                time.parent = NA, time = NA,
                                .rows = nEdges)

  # utility function
  getNodeIndex <- function(i) {
    return(as.integer(nTips + i - 1))
  }

  # Parcourt the table
  counter <- 1
  for (curentHost in 1:nTips) {
    ## All "descendants" of parent
    sub.table <- table.hosts[table.hosts[["inf.by"]] == hosts[curentHost], with = TRUE]

    ## Parent of curent host
    parent <- table.hosts[table.hosts[["indTips"]] == curentHost, "indNodes"]
    t_parent <- table.hosts[table.hosts[["indTips"]] == curentHost, "inf.time"]

    ## Each distinct infection events from the host
    for (child in unique(sub.table[["indNodes"]])) {

      # Informations about the transmission event
      sst <- table.hosts[table.hosts[["indNodes"]] == child, ][1,]
      t_child <- sst[["inf.time"]]

      if (counter == 1) { ## Root
        treeTable[counter, ] <- c(parent = getNodeIndex(child), # Root is its own parent
                                  node = getNodeIndex(child),
                                  branch.length = t_child - t_parent,
                                  label = paste0(sst[["inf.by"]], "_", counter),
                                  host = sst[["inf.by"]],
                                  state = table.hosts[child, "inf.in"],
                                  time.parent = t_parent,
                                  time = t_child)
      } else {
        treeTable[counter, ] <- c(parent = getNodeIndex(parent),
                                  node = getNodeIndex(child),
                                  branch.length = t_child - t_parent,
                                  label = paste0(sst[["inf.by"]], "_", counter),
                                  host = sst[["inf.by"]],
                                  state = table.hosts[child, "inf.in"],
                                  time.parent = t_parent,
                                  time = t_child)
      }

      # actualize
      t_parent <- t_child
      parent <- child
      counter <- counter + 1
    }

    ## Last tip: fictive "dying" host
    t_child <- table.hosts[["out.time"]][curentHost]
    treeTable[counter, ] <- c(parent = getNodeIndex(parent),
                              node = curentHost,
                              branch.length = t_child - t_parent,
                              label = hosts[curentHost],
                              host = hosts[curentHost],
                              state = table.hosts[curentHost, "current.in"],
                              time.parent = t_parent,
                              time = t_child)
    counter <- counter + 1

  }
  # Get correct object
  class(treeTable) <- c("tbl_tree", class(treeTable))
  resTree <- tidytree::as.treedata(treeTable)
  resTree@phylo$root.edge <- treeTable[[1, "branch.length"]]
  return(resTree)
}

## Get infos
get_node <- function(tree, host, time) {
  # Host
  node_bool_host <- (tree@data$host == host)
  if (sum(node_bool_host) < 1) stop(paste0("There are no node with host ", host, " in the tree."))
  # Time
  node_bool <- node_bool_host & (tree@data$time > time) & (tree@data$time.parent <= time)
  if (sum(node_bool) < 1) {
    # Tip case
    node_bool <- node_bool_host & (tree@data$time >= time) & (tree@data$time.parent <= time)
    if (sum(node_bool) < 1) {
      stop(paste0("Host ", host, " is not alive at time ", time, "."))
    }
    if (sum(node_bool) > 1) {
      return(min(tree@data$node[node_bool])) # get the tip if ambiguity
    }
  }
  # Return node number
  return(tree@data$node[node_bool])
}

get_position <- function(tree, node, time) {
  return(tree@data[tree@data$node == node, "time", drop = TRUE] - time)
}

get_state <- function(table.state, host, time, total.time) {
  if (time > total.time) stop(paste0("Time ", time, " is larger than total time ", total.time, " for the epidemic."))
  table.state$time.to[is.na(table.state$time.to)] <- Inf
  node_bool_host <- (table.state$hosts.ID == host)
  if (sum(node_bool_host) < 1) stop(paste0("There are no host named ", host, " in the chain."))
  # Time
  node_bool <- node_bool_host & (table.state$time.to > time) & (table.state$time.from <= time)
  if (sum(node_bool) < 1) {
    stop(paste0("Host ", host, " is not alive at time ", time, "."))
  }
  # Return node number
  return(table.state$state[node_bool])
}

## Add one tip
add_node_tip <- function(tree, node, label, time, host, state) {
  # parameters
  oldTree <- tree@phylo
  oldData <- tidytree::as_tibble(tree)
  # new tip
  tip <- list(edge = matrix(c(2,1),1,2),
              tip.label = label,
              edge.length = 0.0,
              Nnode = 1,
              node.label = paste0("node_", label))
  class(tip) <- "phylo"
  df <- tidytree::tibble(label = c(label, paste0("node_", label)),
                         host = c(host, host),
                         state = c(state, state),
                         time.parent = c(time, oldData[oldData$node == node, "time.parent", drop = TRUE]),
                         time = c(time, time))
  # bind ape objects
  newTree <- ape::bind.tree(oldTree, tip, node, get_position(tree, node, time))
  newTree$node.label[is.na(newTree$node.label)] <- tip$node.label
  newTreeTibble <- tidytree::as_tibble(newTree)
  # bind data objects
  oldData <- oldData[-c(1, 2, 3)]
  newData <- rbind(oldData, df)
  # Create new treedata
  newTibbleTree <- tidytree::full_join(newTreeTibble, newData, by = "label")
  # New Object
  resTree <- tidytree::as.treedata(newTibbleTree)
  resTree@phylo$root.edge <- newTree$root.edge
  return(resTree)
}


## Sample the tree
sampleTransmissionTree <- function(nosoiInf, tree, samples) {
  ## Add all tips
  for (i in 1:nrow(samples)) {
    node <- get_node(tree, samples$hosts[i], samples$times[i])
    state <- get_state(nosoiInf$host.info.A$table.state, samples$hosts[i], samples$times[i], nosoiInf$total.time)

    tree <- add_node_tip(tree,
                         node,
                         samples$labels[i],
                         samples$times[i],
                         samples$hosts[i],
                         state)
  }

  resTree <- treeio::drop.tip(tree, tree@phylo$tip.label[-match(samples$labels, tree@phylo$tip.label)])
  resTree@phylo$root.edge <- tree@phylo$root.edge
  return(resTree)
}

# ## Trim all the tips below a given node
# trim_tree <- function(tree, node, label, relative_time) {
#   # Clade under the node
#   sub_tree <- tree_subset(tree = tree, node = node, levels_back = 0)
#   tips <- get.tree(sub_tree)$tip.label
#   # Add node as tip
#   add_tip(tree)
#   # Trim the whole clade
#   ttree <- drop.tip.treedata(tree, tips, subtree = TRUE)
#   # Label the new tip
#   new_tip_int <- grep("\\[.*\\]", get.tree(ttree)$tip.label)
#   get.tree(ttree)$tip.label[new_tip_int] <- label
#   # Cut to right time
#   new_tip_edge <- get.tree(ttree)$edge[, 2] == new_tip_int
#   get.tree(ttree)$edge.length[new_tip_edge] <- relative_time * get.tree(ttree)$edge.length[new_tip_edge]
#   return(ttree)
# }
