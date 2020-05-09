context("Testing the transmission tree functions")

test_that("Single, discrete", {

  skip_if_not_installed("ape")
  skip_if_not_installed("treeio")
  skip_if_not_installed("tidytree")

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}
  p_Move_fct  <- function(t){return(0.1)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(805)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=20,
                          max.infected=100,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  thostTable <- getTableHosts(test.nosoiA)
  tstateTable <- getTableState(test.nosoiA)
  tot_time <- test.nosoiA$total.time

  ## Full transmission tree
  ttreedata <- getTransmissionTree(test.nosoiA)
  ttree <- ttreedata@phylo
  tdata <- ttreedata@data
  # total time
  expect_equivalent(max(diag(ape::vcv(ttree))) + ttree$root.edge, test.nosoiA$total.time)
  # hosts
  nHosts <- nrow(thostTable)
  for (i in 2:nrow(thostTable)) {
    expect_equivalent(subset(tdata, node == nHosts + thostTable$indNodes[i] - 1)$host,
                      thostTable$inf.by[i])
  }
  # number of descendants
  for (i in 2:nrow(thostTable)) {
    expect_equivalent(length(tidytree::child(ttree, nHosts + thostTable$indNodes[i] - 1)),
                      sum(thostTable$indNodes == thostTable$indNodes[i]) + 1)
  }

  ## Extracting functions
  # Errors
  expect_error(get_node(tdata, "H-123", 2.2), "There are no node with host H-123 in the tree.")
  expect_error(get_node(tdata, "H-7", 2.2), "Host H-7 is not alive at time 2.2.")
  expect_error(get_node(tdata, "H-7", 11.3), "Host H-7 is not alive at time 11.3.")
  #
  expect_error(get_state(tstateTable, "H-123", 2.2, tot_time), "There are no host named H-123 in the chain.")
  expect_error(get_state(tstateTable, "H-7", 2.2, tot_time), "Host H-7 is not alive at time 2.2.")
  expect_error(get_state(tstateTable, "H-7", 11.3, tot_time), "Host H-7 is not alive at time 11.3.")
  expect_error(get_state(tstateTable, "H-7", 16.1, tot_time), "Time 16.1 is larger than total time 16 for the epidemic.")
  # Tip
  expect_equivalent(get_node(tdata, "H-7", 7), 7)
  expect_equivalent(get_node(tdata, "H-7", 10), 7)
  expect_equivalent(get_node(tdata, "H-4", tot_time), 4)
  #
  expect_equivalent(get_position(tdata, 7, 7), 3)
  expect_equivalent(get_position(tdata, 7, 8.5), 1.5)
  expect_equivalent(get_position(tdata, 4, 15.1), 0.9)
  #
  expect_equivalent(get_state(tstateTable, "H-7", 7, tot_time), "A")
  expect_equivalent(get_state(tstateTable, "H-7", 9, tot_time), "C")
  expect_equivalent(get_state(tstateTable, "H-7", 10, tot_time), "C")
  expect_equivalent(get_state(tstateTable, "H-4", tot_time, tot_time), "C")
  # Node
  expect_equivalent(get_node(tdata, "H-1", 6.9), 126)
  expect_equivalent(get_node(tdata, "H-1", 7), 127)
  #
  expect_equivalent(get_position(tdata, 126, 6.5), 0.5)
  expect_equivalent(get_position(tdata, 127, 7.5), 0.5)

  ## Tip states
  for (i in 1:nrow(thostTable)) {
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      thostTable$current.in[i])
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      tdata$state[i])
  }

  ## Node states
  for (i in 1:nrow(thostTable)) {
    for (t in thostTable$inf.time[i]:thostTable$out.time[i]) {
      nn <- get_node(tdata, thostTable$hosts.ID[i], t)
      tt <- tdata$time[tdata$node == nn]
      expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], tt, tot_time),
                        tdata$state[tdata$node == nn])
    }
  }

  ## Full extraction
  hID <- c("H-1", "H-7", "H-15", "H-100")
  samples <- data.table(hosts = hID,
                        times = c(5.2, 9.3, 10.2, 16),
                        labels = paste0(hID, "-s"))

  sampledTree <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  # plot(sampledTree@phylo)

  sttree <- sampledTree@phylo
  stdata <- sampledTree@data
  # total time
  expect_equivalent(max(diag(ape::vcv(sttree))) + sttree$root.edge, test.nosoiA$total.time)
  # DO MORE TESTS

  ## Sampling from the deads
  sampledDeadTree <- sampleTransmissionTreeFromExiting(ttreedata, hID)
  # plot(sampledDeadTree@phylo)
  sampledDeadTreeData <- tidytree::as_tibble(sampledDeadTree)
  sampledDeadTreeData[1:length(hID), ] <- sampledDeadTreeData[match(hID, sampledDeadTreeData$label), ]

  # Check that the two methods give the same results
  samples <- data.table(hosts = hID,
                        times = thostTable[match(hID, thostTable$hosts.ID), "out.time"],
                        labels = paste0(hID, "-s"))
  sampledTreeDeadBis <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  sampledTreeDeadBisData <- tidytree::as_tibble(sampledTreeDeadBis)

  expect_equivalent(sampledDeadTreeData$state, sampledTreeDeadBisData$state)
  expect_equivalent(sampledDeadTreeData$host, sampledTreeDeadBisData$host)
  expect_equivalent(sampledDeadTreeData$time, sampledTreeDeadBisData$time)
  expect_equivalent(sampledDeadTreeData$time.parent, sampledTreeDeadBisData$time.parent)
})

test_that("Single, continuous", {

  skip_if_not_installed("ape")
  skip_if_not_installed("treeio")
  skip_if_not_installed("tidytree")
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
  # plot(test.raster)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  sdMove_fct = function(t,current.env.value){return(100/(current.env.value+1))}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  start.pos <- c(0,0)

  test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                          length=20,
                          max.infected=500,
                          init.individuals=1,
                          init.structure=start.pos,
                          structure.raster=test.raster,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.sdMove=TRUE,
                          sdMove=sdMove_fct,
                          param.sdMove=NA,
                          attracted.by.raster=TRUE,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA)

  thostTable <- getTableHosts(test.nosoiA)
  tstateTable <- getTableState(test.nosoiA)
  tot_time <- test.nosoiA$total.time

  ## Full transmission tree
  ttreedata <- getTransmissionTree(test.nosoiA)
  ttree <- ttreedata@phylo
  tdata <- ttreedata@data
  # total time
  expect_equivalent(max(diag(ape::vcv(ttree))) + ttree$root.edge, test.nosoiA$total.time)
  # hosts
  nHosts <- nrow(thostTable)
  for (i in 2:nrow(thostTable)) {
    expect_equivalent(subset(tdata, node == nHosts + thostTable$indNodes[i] - 1)$host,
                      thostTable$inf.by[i])
  }
  # number of descendants
  for (i in 2:nrow(thostTable)) {
    expect_equivalent(length(tidytree::child(ttree, nHosts + thostTable$indNodes[i] - 1)),
                      sum(thostTable$indNodes == thostTable$indNodes[i]) + 1)
  }

  ## Extracting functions
  expect_equivalent(get_state(tstateTable, "H-5", 7, tot_time),
                    c(state.x = -1.236314, state.y = -1.544381),
                    tolerance = 1.0e-6)
  expect_equivalent(get_state(tstateTable, "H-5", 15, tot_time),
                    c(state.x = -1.615572, state.y = -1.428086),
                    tolerance = 1.0e-6)

  ## Tip states
  for (i in 1:nrow(thostTable)) {
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      c(state.x = thostTable$current.in.x[i], state.y = thostTable$current.in.y[i]),
                      tolerance = 1.0e-6)
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      c(state.x = tdata$state.x[which(tdata$host[1:length(ttree$tip.label)] == thostTable$hosts.ID[i])],
                        state.y = tdata$state.y[which(tdata$host[1:length(ttree$tip.label)] == thostTable$hosts.ID[i])]),
                      tolerance = 1.0e-6)
  }
  ## Node states
  for (i in 1:nrow(thostTable)) {
    for (t in thostTable$inf.time[i]:thostTable$out.time[i]) {
      nn <- get_node(tdata, thostTable$hosts.ID[i], t)
      tt <- tdata$time[tdata$node == nn]
      expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], tt, tot_time),
                        c(state.x = tdata$state.x[tdata$node == nn],
                          state.y = tdata$state.y[tdata$node == nn]),
                        tolerance = 1.0e-6)
    }
  }


  ## Full extraction
  hID <- c("H-1", "H-3", "H-85", "H-5")
  samples <- data.table(hosts = hID,
                        times = c(5.2, 9.3, 20, 10.2),
                        labels = paste0(hID, "-s"))

  sampledTree <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  # plot(sampledTree@phylo)

  sttree <- sampledTree@phylo
  stdata <- sampledTree@data
  # total time
  expect_equivalent(max(diag(ape::vcv(sttree))) + sttree$root.edge, test.nosoiA$total.time)
  # DO MORE TESTS

  ## Sampling from the deads
  sampledDeadTree <- sampleTransmissionTreeFromExiting(ttreedata, hID)
  # plot(sampledDeadTree@phylo)
  sampledDeadTreeData <- tidytree::as_tibble(sampledDeadTree)
  sampledDeadTreeData[1:length(hID), ] <- sampledDeadTreeData[match(hID, sampledDeadTreeData$label), ]

  # Check that the two methods give the same results
  samples <- data.table(hosts = hID,
                        times = thostTable[match(hID, thostTable$hosts.ID), "out.time"],
                        labels = paste0(hID, "-s"))
  sampledTreeDeadBis <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  sampledTreeDeadBisData <- tidytree::as_tibble(sampledTreeDeadBis)

  expect_equivalent(sampledDeadTreeData$state.x, sampledTreeDeadBisData$state.x)
  expect_equivalent(sampledDeadTreeData$state.y, sampledTreeDeadBisData$state.y)
  expect_equivalent(sampledDeadTreeData$host, sampledTreeDeadBisData$host)
  expect_equivalent(sampledDeadTreeData$time, sampledTreeDeadBisData$time)
  expect_equivalent(sampledDeadTreeData$time.parent, sampledTreeDeadBisData$time.parent)

})

test_that("Dual, discrete", {

  skip_if_not_installed("ape")
  skip_if_not_installed("treeio")
  skip_if_not_installed("tidytree")

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  p_Exit_fct  <- function(t,current.in){
    if(current.in=="A"){return(0)}
    if(current.in=="B"){return(0.5)}
    if(current.in=="C"){return(1)}}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(805)
  test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                          length.sim=20,
                          max.infected.A=1000,
                          max.infected.B=1000,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A="A",
                          init.structure.B=NA,
                          structure.matrix.A=transition.matrix,
                          structure.matrix.B=transition.matrix,

                          pExit.A=p_Exit_fct,
                          param.pExit.A=NA,
                          timeDep.pExit.A=FALSE,
                          diff.pExit.A=TRUE,
                          pMove.A=p_Move_fct,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          nContact.A=time_contact,
                          param.nContact.A=NA,
                          timeDep.nContact.A=FALSE,
                          diff.nContact.A=FALSE,
                          pTrans.A=proba,
                          param.pTrans.A=list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          timeDep.pTrans.A=FALSE,
                          diff.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B=p_Exit_fct,
                          param.pExit.B=NA,
                          timeDep.pExit.B=FALSE,
                          diff.pExit.B=TRUE,
                          pMove.B=p_Move_fct,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          nContact.B=time_contact,
                          param.nContact.B=NA,
                          timeDep.nContact.B=FALSE,
                          diff.nContact.B=FALSE,
                          pTrans.B=proba,
                          param.pTrans.B=list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          timeDep.pTrans.B=FALSE,
                          diff.pTrans.B=FALSE,
                          prefix.host.B="V")

  thostTable <- nosoi:::merge_host_tables(test.nosoiA)
  tstateTable <- nosoi:::merge_state_tables(test.nosoiA)
  tot_time <- test.nosoiA$total.time

  thostTable$out.time[is.na(thostTable$out.time)] <- tot_time

  ## Full transmission tree
  ttreedata <- getTransmissionTree(test.nosoiA)
  ttree <- ttreedata@phylo
  tdata <- ttreedata@data
  # total time
  expect_equivalent(max(diag(ape::vcv(ttree))) + ttree$root.edge, test.nosoiA$total.time)

  ## Extracting functions
  # Errors
  expect_error(get_node(tdata, "H-123", 2.2), "There are no node with host H-123 in the tree.")
  expect_error(get_node(tdata, "H-11", 2.2), "Host H-11 is not alive at time 2.2.")
  expect_error(get_node(tdata, "H-11", 14.2), "Host H-11 is not alive at time 14.2.")
  #
  expect_error(get_state(tstateTable, "H-123", 2.2, tot_time), "There are no host named H-123 in the chain.")
  expect_error(get_state(tstateTable, "H-11", 2.2, tot_time), "Host H-11 is not alive at time 2.2.")
  expect_error(get_state(tstateTable, "H-11", 14.2, tot_time), "Host H-11 is not alive at time 14.2.")
  expect_error(get_state(tstateTable, "H-7", 20.1, tot_time), "Time 20.1 is larger than total time 20 for the epidemic.")
  # Tip
  expect_equivalent(get_node(tdata, "H-11", 11), 18)
  expect_equivalent(get_node(tdata, "H-11", 12), 18)
  expect_equivalent(get_node(tdata, "V-66", tot_time), 141)
  #
  expect_equivalent(get_position(tdata, 18, 11), 3)
  expect_equivalent(get_position(tdata, 18, 12.5), 1.5)
  expect_equivalent(get_position(tdata, 141, 15.1), 4.9)
  #
  expect_equivalent(get_state(tstateTable, "H-11", 11, tot_time), "A")
  expect_equivalent(get_state(tstateTable, "H-11", 13, tot_time), "B")
  expect_equivalent(get_state(tstateTable, "H-11", 14, tot_time), "B")
  expect_equivalent(get_state(tstateTable, "V-66", tot_time, tot_time), "A")
  # Node
  expect_equivalent(get_node(tdata, "V-1", 10.9), 150)
  expect_equivalent(get_node(tdata, "V-1", 11), 152)
  #
  expect_equivalent(get_position(tdata, 150, 10.5), 0.5)
  expect_equivalent(get_position(tdata, 152, 11.5), 0.5)

  ## Tip states
  for (i in 1:nrow(thostTable)) {
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      thostTable$current.in[i])
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      tdata$state[which(tdata$host[1:length(ttree$tip.label)] == thostTable$hosts.ID[i])])
  }
  ## Node states
  for (i in 1:nrow(thostTable)) {
    for (t in thostTable$inf.time[i]:thostTable$out.time[i]) {
      nn <- get_node(tdata, thostTable$hosts.ID[i], t)
      tt <- tdata$time[tdata$node == nn]
      expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], tt, tot_time),
                        tdata$state[tdata$node == nn])
    }
  }

  ## Full extraction
  hID <- c("H-1", "H-11", "V-7", "V-66")
  samples <- data.table(hosts = hID,
                        times = c(5.2, 12.6, 6.2, 20),
                        labels = paste0(hID, "-s"))

  sampledTree <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  # plot(sampledTree@phylo)

  sttree <- sampledTree@phylo
  stdata <- sampledTree@data
  # total time
  expect_equivalent(max(diag(ape::vcv(sttree))) + sttree$root.edge, test.nosoiA$total.time)
  # DO MORE TESTS

  ## Sampling from the deads
  sampledDeadTree <- sampleTransmissionTreeFromExiting(ttreedata, hID)
  # plot(sampledDeadTree@phylo)
  sampledDeadTreeData <- tidytree::as_tibble(sampledDeadTree)
  sampledDeadTreeData[1:length(hID), ] <- sampledDeadTreeData[match(hID, sampledDeadTreeData$label), ]

  # Check that the two methods give the same results
  samples <- data.table(hosts = hID,
                        times = thostTable[match(hID, thostTable$hosts.ID), "out.time"],
                        labels = paste0(hID, "-s"))
  samples[is.na(samples$times),"times.out.time"] <- tot_time
  sampledTreeDeadBis <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  sampledTreeDeadBisData <- tidytree::as_tibble(sampledTreeDeadBis)

  expect_equivalent(sampledDeadTreeData$state, sampledTreeDeadBisData$state)
  expect_equivalent(sampledDeadTreeData$host, sampledTreeDeadBisData$host)
  expect_equivalent(sampledDeadTreeData$time, sampledTreeDeadBisData$time)
  expect_equivalent(sampledDeadTreeData$time.parent, sampledTreeDeadBisData$time.parent)
})

test_that("Dual, continuous", {

  skip_if_not_installed("ape")
  skip_if_not_installed("treeio")
  skip_if_not_installed("tidytree")
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  sdMove_fct = function(t,current.env.value){return(100/(current.env.value+1))}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  start.pos <- c(0,0)


  set.seed(805)
  test.nosoiA <- nosoiSim(type="dual", popStructure="continuous",
                          length.sim=200,
                          max.infected.A=500,
                          max.infected.B=500,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A=start.pos,
                          init.structure.B=NA,
                          structure.raster.A=test.raster,
                          structure.raster.B=test.raster,

                          pExit.A=p_Exit_fct,
                          param.pExit.A=NA,
                          timeDep.pExit.A=FALSE,
                          diff.pExit.A=FALSE,
                          pMove.A=p_Move_fct,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          diff.sdMove.A=TRUE,
                          sdMove.A=sdMove_fct,
                          param.sdMove.A=NA,
                          attracted.by.raster.A=TRUE,
                          nContact.A=time_contact,
                          param.nContact.A=NA,
                          timeDep.nContact.A=FALSE,
                          diff.nContact.A=FALSE,
                          pTrans.A=proba,
                          param.pTrans.A=list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          timeDep.pTrans.A=FALSE,
                          diff.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B=p_Exit_fct,
                          param.pExit.B=NA,
                          timeDep.pExit.B=FALSE,
                          diff.pExit.B=FALSE,
                          pMove.B=p_Move_fct,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          diff.sdMove.B=TRUE,
                          sdMove.B=sdMove_fct,
                          param.sdMove.B=NA,
                          attracted.by.raster.B=TRUE,
                          nContact.B=time_contact,
                          param.nContact.B=NA,
                          timeDep.nContact.B=FALSE,
                          diff.nContact.B=FALSE,
                          pTrans.B=proba,
                          param.pTrans.B=list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          timeDep.pTrans.B=FALSE,
                          diff.pTrans.B=FALSE,
                          prefix.host.B="V")


  thostTable <- nosoi:::merge_host_tables(test.nosoiA)
  tstateTable <- nosoi:::merge_state_tables(test.nosoiA)
  tot_time <- test.nosoiA$total.time

  thostTable$out.time[is.na(thostTable$out.time)] <- tot_time

  ## Full transmission tree
  ttreedata <- getTransmissionTree(test.nosoiA)
  ttree <- ttreedata@phylo
  tdata <- ttreedata@data
  # total time
  expect_equivalent(max(diag(ape::vcv(ttree))) + ttree$root.edge, test.nosoiA$total.time)

  ## Extracting functions
  expect_equivalent(get_state(tstateTable, "V-98", 19, tot_time),
                    c(state.x = -1.485684, state.y = -2.014937),
                    tolerance = 1.0e-6)
  expect_equivalent(get_state(tstateTable, "V-98", 21, tot_time),
                    c(state.x = -4.174858, state.y = -1.527998),
                    tolerance = 1.0e-6)

  ## Tip states
  for (i in 1:nrow(thostTable)) {
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      c(state.x = thostTable$current.in.x[i], state.y = thostTable$current.in.y[i]),
                      tolerance = 1.0e-6)
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      c(state.x = tdata$state.x[which(tdata$host[1:length(ttree$tip.label)] == thostTable$hosts.ID[i])],
                        state.y = tdata$state.y[which(tdata$host[1:length(ttree$tip.label)] == thostTable$hosts.ID[i])]),
                      tolerance = 1.0e-6)
  }
  ## Node states
  for (i in 1:nrow(thostTable)) {
    for (t in thostTable$inf.time[i]:thostTable$out.time[i]) {
      nn <- get_node(tdata, thostTable$hosts.ID[i], t)
      tt <- tdata$time[tdata$node == nn]
      expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], tt, tot_time),
                        c(state.x = tdata$state.x[tdata$node == nn],
                          state.y = tdata$state.y[tdata$node == nn]),
                        tolerance = 1.0e-6)
    }
  }

  ## Full extraction
  hID <- c("H-1", "H-3", "H-102", "V-98")
  samples <- data.table(hosts = hID,
                        times = c(5.2, 9.3, 22, 19.3),
                        labels = paste0(hID, "-s"))

  sampledTree <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  # plot(sampledTree@phylo)

  sttree <- sampledTree@phylo
  stdata <- sampledTree@data
  # total time
  expect_equivalent(max(diag(ape::vcv(sttree))) + sttree$root.edge, test.nosoiA$total.time)
  # DO MORE TESTS

  ## Sampling from the deads
  sampledDeadTree <- sampleTransmissionTreeFromExiting(ttreedata, hID)
  # plot(sampledDeadTree@phylo)
  sampledDeadTreeData <- tidytree::as_tibble(sampledDeadTree)
  sampledDeadTreeData[1:length(hID), ] <- sampledDeadTreeData[match(hID, sampledDeadTreeData$label), ]

  # Check that the two methods give the same results
  samples <- data.table(hosts = hID,
                        times = thostTable[match(hID, thostTable$hosts.ID), "out.time"],
                        labels = paste0(hID, "-s"))
  samples[is.na(samples$times),"times.out.time"] <- tot_time
  sampledTreeDeadBis <- sampleTransmissionTree(test.nosoiA, ttreedata, samples)
  sampledTreeDeadBisData <- tidytree::as_tibble(sampledTreeDeadBis)

  expect_equivalent(sampledDeadTreeData$state.x, sampledTreeDeadBisData$state.x)
  expect_equivalent(sampledDeadTreeData$state.y, sampledTreeDeadBisData$state.y)
  expect_equivalent(sampledDeadTreeData$host, sampledTreeDeadBisData$host)
  expect_equivalent(sampledDeadTreeData$time, sampledTreeDeadBisData$time)
  expect_equivalent(sampledDeadTreeData$time.parent, sampledTreeDeadBisData$time.parent)

})

test_that("Single, discrete, Sink", {

  skip_if_not_installed("ape")
  skip_if_not_installed("treeio")
  skip_if_not_installed("tidytree")

  traits <- data.frame(location=rbind('A','B','C', 'D'))
  Q <-list()
  for (column in 1:ncol(traits)) {
    suppressWarnings({ ## We know it fills diagonals with NAs
      Q[[column]] <- diag(unique(traits[,column]), nrow = length(unique(traits[,column])))
    })
    #    diag(Q[[column]]) = 1-nrow(Q[[column]])
    diag(Q[[column]]) = 0
    Q[[column]][lower.tri(Q[[column]])] <- 1/(nrow(Q[[column]])-1)
    Q[[column]][upper.tri(Q[[column]])] <- 1/(nrow(Q[[column]])-1)
    colnames(Q[[column]]) <- rownames(Q[[column]]) <- unique(traits[,column])
  }
  names(Q) <- colnames(traits[-1])

  # #pExit daily probability for a host to leave the simulation (either cured, died, etc.).
  p_Exit_fct  <- function(t, t_sero){
    if(t <= t_sero){p=0}
    if(t > t_sero){p=0.07} #inverse of duration of infection (~14 days)
    return(p)
  }
  t_sero_fct <- function(x){rnorm(x,mean = 14,sd=1)} #Approximately 14 days

  #pMove probability (per unit of time) for a host to do move, i.e. to leave its current state (for example, leaving state “A”). It should not be confused with the probabilities extracted from the structure.matrix, which represent the probability to go to a specific location once a movement is ongoing (for example, going to “B” or “C” while coming from “A”).
  p_Move_fct  <- function(t, prestime, current.in){
    if(current.in=="A" & prestime <= 14){return(0)}
    if(current.in=="A" & prestime > 14){return(0.015)}
    #    return(0.9*1*exp(-0.25*prestime)/((1-0.9)+0.9*exp(-0.25*prestime)))}
    if(current.in=="B"){return(0)}
    if(current.in=="C"){return(0)}
    if(current.in=="D"){return(0)}
  }
  #  t_clust_fct <- function(x){rnorm(x,mean = 3,sd=1)}

  #  nContact is the number (expressed as a positive integer) of potentially infectious contacts an infected hosts can encounter per unit of time.
  n_contact_fct = function(t, current.in) {
    if(current.in=="A"){return(1)}
    if(current.in=="B"){return(3)}
    if(current.in=="C"){return(3)}
    if(current.in=="D"){return(3)}
  }

  #pTrans represents the probability of transmission over time (when a contact occurs).
  # in the form of a threshold function: before a certain amount of time since initial infection, the host does not transmit (incubation time, which we call t_incub), and after that time, it will transmit with a certain (constant) probability (which we call p_max). This function is dependent of the time since the host’s infection t.
  p_Trans_fct <- function(t, prestime, current.in, t_incub, p_max){
    if(t < t_incub){p=0}
    if(t >= t_incub & current.in=="A"){p=p_max}
    if(t >= t_incub & current.in=="B"){p=0.2}
    if(t >= t_incub & current.in=="C"){p=0.8*1*exp(-0.1*prestime)/((1-0.8)+0.8*exp(-0.1*prestime))} ##0.9 is the initial, 1 is the carrying capacity, -0.5 is the growth rate
    if(t >= t_incub & current.in=="D"){p=0.2*1*exp(0.0025*prestime)/((1-0.2)+0.2*exp(0.0025*prestime))}
    return(p)
  }

  t_incub_fct <- function(x){rnorm(x,mean = 4,sd=1)} #Approximately 4.2 days
  p_max_fct <- function(x){rbeta(x,shape1 = 1,shape2=5)}# Mean of 0.14

  # Starting the simulation ------------------------------------

  set.seed(858)

  SimulationSingle <- nosoiSim(type="single", # Number of hosts
                               popStructure="discrete", #discrete or continuous
                               structure.matrix = Q[[1]], # prob matrix defined above (row sums to 1, with diags as zero)
                               length.sim = 400, # Max number of time units (can be days, months, weeks, etc.)
                               max.infected = 40, #maximum number of individuals that can be infected during the simulation.
                               init.individuals = 1, #number of individuals (an integer above 1) that will start a transmission chain. Keep in mind that you will have as many transmission chains as initial individuals, which is equivalent as launching a number of independent nosoi simulations.
                               init.structure = "A",

                               pExit = p_Exit_fct,
                               param.pExit=list(t_sero=t_sero_fct),
                               timeDep.pExit=FALSE,
                               diff.pExit=FALSE,

                               pMove = p_Move_fct,
                               param.pMove=NA,
                               timeDep.pMove=TRUE,
                               diff.pMove=TRUE,

                               nContact=n_contact_fct,
                               param.nContact = NA,
                               timeDep.nContact=FALSE,
                               diff.nContact=TRUE,

                               # nContact=time_contact,
                               # hostCount.nContact=TRUE,
                               # param.nContact = NA,
                               # timeDep.nContact=FALSE,
                               # diff.nContact=TRUE,

                               pTrans = p_Trans_fct,
                               param.pTrans = list(p_max=p_max_fct, t_incub=t_incub_fct),
                               timeDep.pTrans=TRUE,
                               diff.pTrans=TRUE,

                               prefix.host="C",
                               print.progress=TRUE,
                               print.step=10)

  thostTable <- getTableHosts(SimulationSingle)
  tstateTable <- getTableState(SimulationSingle)
  tot_time <- SimulationSingle$total.time

  ## Full transmission tree
  ttreedata <- getTransmissionTree(SimulationSingle)
  ttree <- ttreedata@phylo
  tdata <- ttreedata@data
  # total time
  expect_equivalent(max(diag(ape::vcv(ttree))) + ttree$root.edge, SimulationSingle$total.time)
  # hosts
  nHosts <- nrow(thostTable)
  for (i in 2:nrow(thostTable)) {
    expect_equivalent(subset(tdata, node == nHosts + thostTable$indNodes[i] - 1)$host,
                      thostTable$inf.by[i])
  }
  # number of descendants
  for (i in 2:nrow(thostTable)) {
    expect_equivalent(length(tidytree::child(ttree, nHosts + thostTable$indNodes[i] - 1)),
                      sum(thostTable$indNodes == thostTable$indNodes[i]) + 1)
  }

  ## Tip states
  for (i in 1:nrow(thostTable)) {
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      thostTable$current.in[i])
    expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], thostTable$out.time[i], tot_time),
                      tdata$state[i])
  }

  ## Node states
  for (i in 1:nrow(thostTable)) {
    for (t in thostTable$inf.time[i]:thostTable$out.time[i]) {
      nn <- get_node(tdata, thostTable$hosts.ID[i], t)
      tt <- tdata$time[tdata$node == nn]
      expect_equivalent(get_state(tstateTable, thostTable$hosts.ID[i], tt, tot_time),
                        tdata$state[tdata$node == nn])
    }
  }

})
