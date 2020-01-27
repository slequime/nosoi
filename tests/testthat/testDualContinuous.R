context("Testing dual-host with continuous structure")

test_that("Both hosts move", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  skip_if_not_installed("igraph")
  library(igraph)

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

  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts,test.nosoiA$host.info.B$table.hosts))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 8)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 22)

  expect_equal(test.nosoiA$host.info.A$N.infected, 326)
  expect_equal(test.nosoiA$host.info.B$N.infected, 579)

  expect_equal(test.nosoiA$type, "dual")
  expect_equal(test.nosoiA$host.info.A$popStructure, "continuous")
  expect_equal(test.nosoiA$host.info.B$popStructure, "continuous")

  #Movement
  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "H-1")),2)

})


test_that("One host (A) moves", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  skip_if_not_installed("igraph")
  library(igraph)

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
                          pMove.B=NA,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          diff.sdMove.B=FALSE,
                          sdMove.B=NA,
                          param.sdMove.B=NA,
                          attracted.by.raster.B=FALSE,
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

  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts,test.nosoiA$host.info.B$table.hosts))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 10)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)


  expect_equal(nrow(test.nosoiA$host.info.B$table.hosts), nrow(test.nosoiA$host.info.B$table.state))

  expect_equal(test.nosoiA$total.time, 24)

  expect_equal(test.nosoiA$host.info.A$N.infected, 682)
  expect_equal(test.nosoiA$host.info.B$N.infected, 606)

  expect_equal(test.nosoiA$type, "dual")
  expect_equal(test.nosoiA$host.info.A$popStructure, "continuous")
  expect_equal(test.nosoiA$host.info.B$popStructure, "continuous")

  #Movement

  H1_moves <- subset(full.results.nosoi.state, hosts.ID == "H-1")

  expect_equal(nrow(H1_moves),5)
  expect_equal(H1_moves$current.env.value[1] < H1_moves$current.env.value[5],TRUE)
})

test_that("One host (B) moves", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  skip_if_not_installed("igraph")
  library(igraph)

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

  set.seed(19)
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
                          pMove.A=NA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          diff.sdMove.A=TRUE,
                          sdMove.A=NA,
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
                          attracted.by.raster.B=FALSE,
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

  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts,test.nosoiA$host.info.B$table.hosts))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 10)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)


  expect_equal(nrow(test.nosoiA$host.info.A$table.hosts), nrow(test.nosoiA$host.info.A$table.state))

  expect_equal(test.nosoiA$total.time, 26)

  expect_equal(test.nosoiA$host.info.A$N.infected, 627)
  expect_equal(test.nosoiA$host.info.B$N.infected, 520)

  expect_equal(test.nosoiA$type, "dual")
  expect_equal(test.nosoiA$host.info.A$popStructure, "continuous")
  expect_equal(test.nosoiA$host.info.B$popStructure, "continuous")

  #Movement

  H1_moves <- subset(full.results.nosoi.state, hosts.ID == "V-1")

  expect_equal(nrow(H1_moves),2)
  expect_equal(H1_moves$current.env.value[1] < H1_moves$current.env.value[2],TRUE)
})

test_that("Epidemic dies out", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  skip_if_not_installed("igraph")
  library(igraph)

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

  set.seed(1000)
  test.nosoiA <- nosoiSim(type="dual", popStructure="continuous",
                          length.sim=200,
                          max.infected.A=500,
                          max.infected.B=500,
                          init.individuals.A=0,
                          init.individuals.B=1,
                          init.structure.A=NA,
                          init.structure.B=start.pos,
                          structure.raster.A=test.raster,
                          structure.raster.B=test.raster,

                          pExit.A=p_Exit_fct,
                          param.pExit.A=NA,
                          timeDep.pExit.A=FALSE,
                          diff.pExit.A=FALSE,
                          pMove.A=NA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          diff.sdMove.A=TRUE,
                          sdMove.A=NA,
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
                          attracted.by.raster.B=FALSE,
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

  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts,test.nosoiA$host.info.B$table.hosts))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  expect_equal(nrow(test.nosoiA$host.info.A$table.hosts), nrow(test.nosoiA$host.info.A$table.state))

  expect_equal(test.nosoiA$total.time, 4)

  expect_equal(test.nosoiA$host.info.A$N.infected, 0)

  expect_equal(test.nosoiA$host.info.B$N.infected, 1)

  expect_equal(test.nosoiA$type, "dual")
  expect_equal(test.nosoiA$host.info.A$popStructure, "continuous")
  expect_equal(test.nosoiA$host.info.B$popStructure, "continuous")

  #Movement

  H1_moves <- subset(full.results.nosoi.state, hosts.ID == "V-1")
  expect_equal(nrow(H1_moves),2)
})

test_that("Error if no host move", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  skip_if_not_installed("igraph")
  library(igraph)

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

  set.seed(1000)
  expect_error(
    nosoiSim(type="dual", popStructure="continuous",
             length.sim=200,
             max.infected.A=500,
             max.infected.B=500,
             init.individuals.A=0,
             init.individuals.B=1,
             init.structure.A=NA,
             init.structure.B=start.pos,
             structure.raster.A=test.raster,
             structure.raster.B=test.raster,

             pExit.A=p_Exit_fct,
             param.pExit.A=NA,
             timeDep.pExit.A=FALSE,
             diff.pExit.A=FALSE,
             pMove.A=NA,
             param.pMove.A=NA,
             timeDep.pMove.A=FALSE,
             diff.pMove.A=FALSE,
             diff.sdMove.A=TRUE,
             sdMove.A=NA,
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
             pMove.B=NA,
             param.pMove.B=NA,
             timeDep.pMove.B=FALSE,
             diff.pMove.B=FALSE,
             diff.sdMove.B=TRUE,
             sdMove.B=NA,
             param.sdMove.B=NA,
             attracted.by.raster.B=FALSE,
             nContact.B=time_contact,
             param.nContact.B=NA,
             timeDep.nContact.B=FALSE,
             diff.nContact.B=FALSE,
             pTrans.B=proba,
             param.pTrans.B=list(p_max=p_max_fct,
                                 t_incub=t_incub_fct),
             timeDep.pTrans.B=FALSE,
             diff.pTrans.B=FALSE,
             prefix.host.B="V"),
    "At least one host must move.")
})


test_that("One host (B) moves, host count", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  skip_if_not_installed("igraph")
  library(igraph)

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

  time_contact_A = function(t){round(rnorm(1, 3, 1), 0)}

  time_contact_B <- function(t, current.env.value, host.count.A){

    temp.val = round(((current.env.value-host.count.A)/current.env.value)*rnorm(1, 3, 1), 0)

    if(length(temp.val) == 0 || temp.val <= 0) {
      return(0)
    }
    if(temp.val >= 0) {
      return(temp.val)
    }
  }

  start.pos <- c(0,0)

  set.seed(19)
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
                          pMove.A=NA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          diff.sdMove.A=TRUE,
                          sdMove.A=NA,
                          param.sdMove.A=NA,
                          attracted.by.raster.A=TRUE,
                          nContact.A=time_contact_A,
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
                          attracted.by.raster.B=FALSE,
                          nContact.B=time_contact_B,
                          param.nContact.B=NA,
                          timeDep.nContact.B=FALSE,
                          hostCount.nContact.B=TRUE,
                          diff.nContact.B=TRUE,
                          pTrans.B=proba,
                          param.pTrans.B=list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          timeDep.pTrans.B=FALSE,
                          diff.pTrans.B=FALSE,
                          prefix.host.B="V")

  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts,test.nosoiA$host.info.B$table.hosts))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 9)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)


  expect_equal(nrow(test.nosoiA$host.info.A$table.hosts), nrow(test.nosoiA$host.info.A$table.state))

  expect_equal(test.nosoiA$total.time, 25)

  expect_equal(test.nosoiA$host.info.A$N.infected, 280)
  expect_equal(test.nosoiA$host.info.B$N.infected, 540)

  expect_equal(test.nosoiA$type, "dual")
  expect_equal(test.nosoiA$host.info.A$popStructure, "continuous")
  expect_equal(test.nosoiA$host.info.B$popStructure, "continuous")

  #Movement

  H1_moves <- subset(full.results.nosoi.state, hosts.ID == "V-1")

  expect_equal(nrow(H1_moves),6)
  expect_equal(H1_moves$current.env.value[1] < H1_moves$current.env.value[2],TRUE)
})
