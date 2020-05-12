#' @title Gets the full transmission tree (phylogenetic tree-like) from a \code{nosoi} simulation
#'
#' @description
#'  From a \code{nosoi} simulated epidemics, this function extracts the full transmission tree in a form mimicking a phylogenetic tree.
#'
#' @details
#'  This function uses packages \pkg{tidytree} and \pkg{treeio},
#'  that rely on \code{\link[ape:ape-package]{ape}}.
#'
#' @param nosoiInf an object of class \code{\link{nosoiSim}}
#'
#' @return A tree of class \code{\link[tidytree:treedata-class]{treedata}}, containing a
#' phylogenetic tree based on the transmission chain and the mapped data at all the nodes.
#'
#' @examples
#' \donttest{
#' t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#' p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#' p_Exit_fct  <- function(t){return(0.08)}
#' p_Move_fct  <- function(t){return(0.1)}
#'
#' proba <- function(t,p_max,t_incub){
#'   if(t <= t_incub){p=0}
#'   if(t >= t_incub){p=p_max}
#'   return(p)
#' }
#'
#' time_contact = function(t){round(rnorm(1, 3, 1), 0)}
#'
#' transition.matrix = matrix(c(0, 0.2, 0.4, 0.5, 0, 0.6, 0.5, 0.8, 0),
#'                            nrow = 3, ncol = 3,
#'                            dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#'
#' set.seed(805)
#' test.nosoi <- nosoiSim(type="single", popStructure="discrete",
#'                         length=20,
#'                         max.infected=100,
#'                         init.individuals=1,
#'                         init.structure="A",
#'                         structure.matrix=transition.matrix,
#'                         pMove=p_Move_fct,
#'                         param.pMove=NA,
#'                         nContact=time_contact,
#'                         param.nContact=NA,
#'                         pTrans = proba,
#'                         param.pTrans = list(p_max=p_max_fct,
#'                                             t_incub=t_incub_fct),
#'                         pExit=p_Exit_fct,
#'                         param.pExit=NA
#' )
#'
#' ## Make sure all needed packages are here
#' if (requireNamespace("ape", quietly = TRUE)
#'     || requireNamespace("tidytree", quietly = TRUE)
#'     || requireNamespace("treeio", quietly = TRUE)) {
#'   library(ape)
#'   library(tidytree)
#'   library(treeio)
#'
#'   #' ## Full transmission tree
#'   ttreedata <- getTransmissionTree(test.nosoi)
#'   plot(ttreedata@phylo)
#'
#'   ## Sampling "non dead" individuals
#'   hID <- c("H-1", "H-7", "H-15", "H-100")
#'   samples <- data.table(hosts = hID,
#'                         times = c(5.2, 9.3, 10.2, 16),
#'                         labels = paste0(hID, "-s"))
#'
#'   sampledTree <- sampleTransmissionTree(test.nosoi, ttreedata, samples)
#'   plot(sampledTree@phylo)
#'
#'   ## Sampling "dead" individuals
#'   sampledDeadTree <- sampleTransmissionTreeFromExiting(ttreedata, hID)
#'   plot(sampledDeadTree@phylo)
#'   }
#' }
#'
#' @seealso For exporting the annotated tree to other software packages, see functions
#' in \pkg{treeio} (e.g. \code{\link[treeio:write.beast]{write.beast}}).
#'
#' To sub-sample this tree, see functions \code{\link{sampleTransmissionTree}} and \code{\link{sampleTransmissionTreeFromExiting}}
#'
#' @export getTransmissionTree

getTransmissionTree <- function(nosoiInf) {
  if (!requireNamespace("ape", quietly = TRUE) || !requireNamespace("tidytree", quietly = TRUE) || !requireNamespace("treeio", quietly = TRUE)) {
    stop("Packages 'ape', 'tidytree' and 'treeio' are needed for transmission tree generation.",
         call. = FALSE)
  }
  #To avoid notes (use of dplyr)
  # node <- NULL

  table.hosts <- merge_host_tables(nosoiInf)
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
  treeTable <- tidytree::tibble(parent = NA_integer_, node = NA_integer_, branch.length = NA_real_,
                                label = NA_character_, host = NA_character_,
                                state = NA_character_, state.x = NA_real_, state.y = NA_real_,
                                time.parent = NA_real_, time = NA_real_,
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
                                  state = safeGet(sst, 1, "inf.in"),
                                  state.x = safeGet(sst, 1, "inf.in.x"),
                                  state.y = safeGet(sst, 1, "inf.in.y"),
                                  time.parent = t_parent,
                                  time = t_child)
      } else {
        treeTable[counter, ] <- c(parent = getNodeIndex(parent),
                                  node = getNodeIndex(child),
                                  branch.length = t_child - t_parent,
                                  label = paste0(sst[["inf.by"]], "_", counter),
                                  host = sst[["inf.by"]],
                                  state = safeGet(sst, 1, "inf.in"),
                                  state.x = safeGet(sst, 1, "inf.in.x"),
                                  state.y = safeGet(sst, 1, "inf.in.y"),
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
                              state = safeGet(table.hosts, curentHost, "current.in"),
                              state.x = safeGet(table.hosts, curentHost, "current.in.x"),
                              state.y = safeGet(table.hosts, curentHost, "current.in.y"),
                              time.parent = t_parent,
                              time = t_child)
    counter <- counter + 1

  }
  # Remove NAs
  popStructure <- getHostData(nosoiInf, "popStructure")
  switch(popStructure,
         discrete = treeTable$state.x <- treeTable$state.y <- NULL,
         continuous = treeTable$state <- NULL)

  # Get correct object
  root_length <- treeTable[[1, "branch.length"]]
  # treeTable <- dplyr::arrange(treeTable, node)
  treeTable <- treeTable[order(treeTable[["node"]]), ]
  class(treeTable) <- c("tbl_tree", class(treeTable))
  resTree <- tidytree::as.treedata(treeTable)
  resTree@phylo$root.edge <- root_length
  # resTree@info <- list(pop = pop)
  return(resTree)
}

#' @title Merge Population Data
#'
#' @description
#'  Utility function to get merged hosts data if dual.
#'
#' @param nosoiInf an object of class \code{\link{nosoiSim}}
#'
#' @return the merged `table.hosts` or `table.state` of the two populations.
#'
#' @keywords internal
#'
merge_host_tables <- function(nosoiInf) {
  if (nosoiInf$type == "single") {
    return(nosoiInf$host.info.A$table.hosts) # Just the bare table
  } else if (nosoiInf$type == "dual") {
    # table_A <- getTableHosts(nosoiInf, pop = target_pop)
    # table_B <- getTableHosts(nosoiInf, pop = vector_pop)
    # table_A[table_B, on="inf.by == hosts.ID", nomatch = NULL] # For keeping only one host
    # return(merge(table_A, table_B, by = intersect(colnames(table_A), colnames(table_B)), all = TRUE))
    return(rbindlist(list(nosoiInf$host.info.A$table.hosts, nosoiInf$host.info.B$table.hosts),fill=TRUE))
  }
  stop("Type of transmission should be 'single' or 'dual'-host.")
}


#' @rdname merge_host_tables
merge_state_tables <- function(nosoiInf) {
  if (nosoiInf$type == "single") {
    return(nosoiInf$host.info.A$table.state) # Just the bare table
  } else if (nosoiInf$type == "dual") {
    return(rbindlist(list(nosoiInf$host.info.A$table.state, nosoiInf$host.info.B$table.state)))
  }
  stop("Type of transmission should be 'single' or 'dual'-host.")
}

## Utility functions to get entries in the table, returning NA if does not exist.
is.error <- function(x) inherits(x, "try-error")
safeGet <- function(dt, i, name) {
  res <- try(dt[i, name, with = FALSE], silent = TRUE)
  if (is.error(res)) return(NA);
  return(res)
}

#' @title Get Node
#'
#' @description
#'  Find the node descending from the branch where the host lies at a given
#'  time on the transmission tree.
#'
#' @param tdata a tibble extracted from a \code{treedata} object.
#' @param host host id
#' @param time time at which to sample the host
#'
#' @return number of the node in the transmission tree descending from the
#' host at given time.
#'
#' @seealso \code{\link{get_position}}
#'
#' @keywords internal
##
get_node <- function(tdata, host, time) {
  # Host
  node_bool_host <- (tdata$host == host)
  if (sum(node_bool_host) < 1) stop(paste0("There are no node with host ", host, " in the tree."))
  # Time
  node_bool <- node_bool_host & (tdata$time > time) & (tdata$time.parent <= time)
  if (sum(node_bool) < 1) {
    # Tip case
    node_bool <- node_bool_host & (tdata$time >= time) & (tdata$time.parent <= time)
    if (sum(node_bool) < 1) {
      stop(paste0("Host ", host, " is not alive at time ", time, "."))
    }
    if (sum(node_bool) > 1) {
      return(min(tdata$node[node_bool])) # get the tip if ambiguity
    }
  }
  # Return node number
  return(tdata$node[node_bool])
}

#' @title Get Position on branch
#'
#' @description
#'  Find the position on the branch above the node of the times sample.
#'  Warning: The node needs to be extracted with \code{\link{get_node}}.
#'  Result of this function is to be used in \code{\link[ape]{bind.tree}}
#'
#' @param tdata a tibble extracted from a \code{treedata} object.
#' @param node node below the sampling event
#' @param time time of the sampling event
#'
#' @return time between the node and the sampling event.
#'
#' @seealso \code{\link{get_node}}
#'
#' @keywords internal
##
get_position <- function(tdata, node, time) {
  return(tdata[tdata$node == node, "time", drop = TRUE] - time)
}

#' @title Get State at sampling time
#'
#' @description
#'  Find the state of the host at the sampling time
#'
#' @param table.state data.table of hosts movement, extracted from a \code{nosoi} object
#' @param host ID of the host
#' @param time time of the sampling event
#' @param total.time total time of the epidemics, extracted from the \code{nosoi} object
#'
#' @return state of the host
#'
#' @keywords internal
##
get_state <- function(table.state, host, time, total.time) {
  if (time > total.time) stop(paste0("Time ", time, " is larger than total time ", total.time, " for the epidemic."))
  table.state$time.to[is.na(table.state$time.to)] <- Inf
  node_bool_host <- (table.state$hosts.ID == host)
  if (sum(node_bool_host) < 1) stop(paste0("There are no host named ", host, " in the chain."))
  # Time
  node_bool <- node_bool_host & (table.state$time.to >= time) & (table.state$time.from <= time)
  if (sum(node_bool) < 1) {
    stop(paste0("Host ", host, " is not alive at time ", time, "."))
  }
  # Return node number
  return(c(state = table.state$state[node_bool][length(table.state$state[node_bool])], # get last row if ambiguity
           state.x = table.state$state.x[node_bool][length(table.state$state.x[node_bool])],
           state.y = table.state$state.y[node_bool][length(table.state$state.y[node_bool])]))
}

## Add one tip
#' @title Add one tip
#'
#' @description
#'  Add a tip to the transmission tree corresponding to the sampled individual.
#'
#' @param tree transmission tree of class \code{treedata}, result of \code{\link{getTransmissionTree}}
#' @param host ID of the host
#' @param time time of the sampling event
#' @param label label of the new tip to be added
#' @param state state of the sampled individual
#'
#' @return modified tree, with on extra tip of length zero at the right place of the tree
#'
#' @keywords internal
##
add_node_tip <- function(tree, host, time, label, state) {
  # Tree info
  node <- get_node(tree@data, host, time)
  oldTree <- tree@phylo
  oldData <- tidytree::as_tibble(tree)
  # Check label
  if (label %in% c(oldTree$tip.label, oldTree$node.label)) {
    stop(paste0("Label ", label, " is invalid: please choose a label that is unique."))
  }
  # new tip
  tip <- list(edge = matrix(c(2,1),1,2),
              tip.label = label,
              edge.length = 0.0,
              Nnode = 1,
              node.label = paste0("node_", label))
  class(tip) <- "phylo"
  df <- tidytree::tibble(label = c(label, paste0("node_", label)),
                         host = c(host, host),
                         time.parent = c(oldData[oldData$node == node, "time.parent", drop = TRUE], oldData[oldData$node == node, "time.parent", drop = TRUE]),
                         time = c(time, time))
  if (length(state) == 1) {
    df$state <- c(state, state)
  } else {
    df$state.x <- c(state["state.x"], state["state.x"])
    df$state.y <- c(state["state.y"], state["state.y"])
  }
  # bind ape objects
  position_tip <- get_position(tree@data, node, time)
  if (position_tip == 0 && node <= length(oldTree$tip.label)) {
    ## Special case when trying to add a tip to a tip
    ed <- which(oldTree$edge[, 2] == node)
    ed_l <- oldTree$edge.length[ed]
    node_label <- oldTree$tip.label[node]
    oldTree$edge.length[ed] <- ed_l + 1 ## Add one at the tip
    newTree <- ape::bind.tree(oldTree, tip, node, 1) ## Bind the new tip there
    new_node <- which(newTree$tip.label == node_label)
    newTree$edge.length[which(newTree$edge[, 2] == new_node)] <- ed_l - 1 ## Remove one ot the tip
  } else {
    newTree <- ape::bind.tree(oldTree, tip, node, position_tip)
  }
  newTree$node.label[is.na(newTree$node.label)] <- tip$node.label
  newTreeTibble <- tidytree::as_tibble(newTree)
  # bind data objects
  oldData <- oldData[-c(1, 2, 3)]
  newData <- rbind(df, oldData)
  # Create new treedata
  newTibbleTree <- tidytree::full_join(newTreeTibble, newData, by = "label")
  # New Object
  resTree <- tidytree::as.treedata(newTibbleTree)
  resTree@phylo$root.edge <- newTree$root.edge
  return(resTree)
}

## Vectorized version
draw_one_sample <- function(table.states, total.time, tree, sample) {
  state <- get_state(table.states, sample$hosts, sample$times, total.time)
  return(add_node_tip(tree,
                      sample$hosts,
                      sample$times,
                      sample$labels,
                      state))
}

#' @title Sample the transmission tree (phylogenetic tree-like)
#'
#' @description
#'  Sample a full transmission tree. This function allows for sampling multiple
#'  times on the same lineage. When this happens, the sampled ancestor is
#'  a tip with length zero.
#'
#' @details
#'  The tree needs to be produced by function \code{\link{getTransmissionTree}}
#'  applied on the same \code{nosoiSim} object.
#'
#' @param nosoiInf an object of class \code{\link{nosoiSim}}
#' @param tree a \code{\link[tidytree:treedata]{treedata}} object created by function \code{\link{getTransmissionTree}}
#' @param samples a \code{\link[data.table:data.table-package]{data.table}} object with the following entries:
#' \describe{
#'   \item{hosts}{Host ID of the individuals to be sampled}
#'   \item{times}{Times at which each host is sampled}
#'   \item{labels}{label for the corresponding tip in the tree}
#' }
#'
#' @return A tree of class \code{\link[tidytree:treedata-class]{treedata}}, containing a
#' phylogenetic tree based on the transmission chain and the mapped data at all the nodes.
#'
#' @inherit getTransmissionTree examples
#'
#' @seealso For exporting the annotated tree to other software packages, see functions
#' in \pkg{treeio} (e.g. \code{\link[treeio:write.beast]{write.beast}}).
#'
#' To get the full transmission matrix, see \code{\link{getTransmissionTree}}.
#'
#' For sampling only dead individuals, see \code{\link{sampleTransmissionTreeFromExiting}}.
#'
#' @export sampleTransmissionTree

sampleTransmissionTree <- function(nosoiInf, tree, samples) {
  ## Extract table state
  # pop <- tree@info$pop
  # if (is.null(pop) || !(pop == "A" || pop == "B")){
  #   stop("The tree object has an incorect format. It should be produced by function 'getTransmissionTree'. See documentation.")
  # }
  table.states <- merge_state_tables(nosoiInf)
  tottime <- nosoiInf$total.time
  ## Add all tips
  for (i in 1:nrow(samples)) {
    tree <- draw_one_sample(table.states, tottime, tree, samples[i, ])
  }

  # resTree <- treeio::drop.tip(tree, tree@phylo$tip.label[-match(samples$labels, tree@phylo$tip.label)])
  # resTree@phylo$root.edge <- tree@phylo$root.edge
  resTree <- keep.tip.treedata(tree, samples$labels)
  return(resTree)
}

#' @title Sample the transmission tree (phylogenetic tree-like) among the exited hosts
#'
#' @description
#'  Sample a full transmission tree. This function allows for sampling only exited (i.e. inactive)
#'  individuals (e.g. when the sampling procedure is destructive or cuts the hosts from the population).
#'  Beware because it does not influence the epidemiological process, it only means that the host
#'  has been sampled when exiting the simulation.
#'
#' @details
#'  The tree needs to be produced by function \code{\link{getTransmissionTree}}
#'  applied on the same \code{nosoiSim} object.
#'
#' @param tree a \code{\link[tidytree:treedata-class]{treedata}} object created by function \code{\link{getTransmissionTree}}
#' @param hosts a vector of dead hosts to sample
#'
#' @return A tree of class \code{\link[tidytree:treedata-class]{treedata}}, containing a
#' phylogenetic tree based on the transmission chain and the mapped data at all the nodes.
#'
#' @inherit getTransmissionTree examples
#'
#' @seealso For exporting the annotated tree to other software packages, see functions
#' in \pkg{treeio} (e.g. \code{\link[treeio:write.beast]{write.beast}}).
#'
#' To get the full transmission matrix, see \code{\link{getTransmissionTree}}.
#'
#' For sampling among non-dead individuals, see \code{\link{sampleTransmissionTree}}.
#'
#' @export sampleTransmissionTreeFromExiting

sampleTransmissionTreeFromExiting <- function(tree, hosts) {
  # resTree <- treeio::drop.tip(tree, tree@phylo$tip.label[-match(hosts, tree@phylo$tip.label)])
  # resTree@phylo$root.edge <- tree@phylo$root.edge
  resTree <- keep.tip.treedata(tree, hosts)
  return(resTree)
}

#' @title  Keep tips
#'
#' @description
#'  Keep the tips in the list. See \code{\link[ape:keep.tip]{keep.tip}}
#'
#' @keywords internal

keep.tip.treedata <- function(tree, tip) {
  # Tree info
  oldTree <- tree@phylo
  oldData <- tidytree::as_tibble(tree)
  # keep.tip ape objects
  newTree <- ape::keep.tip(oldTree, tip)
  # mrca
  pos <- oldTree$edge[, 2] == ape::getMRCA(oldTree, tip)
  if (any(pos)) {
    root_length <- oldTree$edge.length[pos]
  } else {
    root_length <- oldTree$root.edge
  }
  # tible
  newTreeTibble <- tidytree::as_tibble(newTree)
  # bind data objects
  oldData <- oldData[-c(1, 2, 3)]
  # Create new treedata
  # newTibbleTree <- dplyr::left_join(newTreeTibble, oldData, by = "label")
  newTibbleTree <- merge(newTreeTibble, oldData, by = "label", all.x = TRUE, sort = FALSE)
  newTibbleTree <- tidytree::as_tibble(newTibbleTree)
  class(newTibbleTree) <- c("tbl_tree", class(newTibbleTree))
  # New Object
  resTree <- tidytree::as.treedata(newTibbleTree)
  resTree@phylo$root.edge <- root_length
  return(resTree)
}
