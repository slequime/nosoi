context("Testing dual-host with continuous structure")

test_that("Both hosts move", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

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
  test.nosoiA <- nosoiSim(type="dual",structure=TRUE, continuous = TRUE,
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
                          prefix.host.B="V",

                          progress.bar=TRUE,
                          print.step=10)

  full.results.nosoi <- rbindlist(list(test.nosoiA$table.hosts_A,test.nosoiA$table.hosts_B))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$table.state_A,test.nosoiA$table.state_B))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 8)

  expect_equal(all(str_detect(test.nosoiA$table.hosts_A$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_A[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 22)

  expect_equal(test.nosoiA$N.infected_A, 326)
  expect_equal(test.nosoiA$N.infected_B, 579)

  expect_equal(test.nosoiA$type, "dualContinuous")

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
test.nosoiA <- nosoiSim(type="dual",structure=TRUE, continuous = TRUE,
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
                        prefix.host.B="V",

                        progress.bar=TRUE,
                        print.step=10)

 full.results.nosoi <- rbindlist(list(test.nosoiA$table.hosts_A,test.nosoiA$table.hosts_B))
 full.results.nosoi.state <- rbindlist(list(test.nosoiA$table.state_A,test.nosoiA$table.state_B))

g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

expect_equal(transitivity(g, type="global"), 0)
expect_equal(clusters(g, "weak")$no, 1)
expect_equal(diameter(g, directed=F, weights=NA), 10)

expect_equal(all(str_detect(test.nosoiA$table.hosts_A$inf.by,"H-") == FALSE),TRUE)
expect_equal(all(str_detect(test.nosoiA$table.hosts_A[-1]$inf.by,"V-") == TRUE),TRUE)
expect_equal(all(str_detect(test.nosoiA$table.hosts_B$inf.by,"V-") == FALSE),TRUE)
expect_equal(all(str_detect(test.nosoiA$table.hosts_B[-1]$inf.by,"H-") == TRUE),TRUE)


expect_equal(nrow(test.nosoiA$table.hosts_B), nrow(test.nosoiA$table.state_B))

expect_equal(test.nosoiA$total.time, 24)

expect_equal(test.nosoiA$N.infected_A, 682)
expect_equal(test.nosoiA$N.infected_B, 606)

expect_equal(test.nosoiA$type, "dualContinuous")

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
  test.nosoiA <- nosoiSim(type="dual",structure=TRUE, continuous = TRUE,
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
                          prefix.host.B="V",

                          progress.bar=TRUE,
                          print.step=10)

  full.results.nosoi <- rbindlist(list(test.nosoiA$table.hosts_A,test.nosoiA$table.hosts_B))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$table.state_A,test.nosoiA$table.state_B))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 10)

  expect_equal(all(str_detect(test.nosoiA$table.hosts_A$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_A[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B[-1]$inf.by,"H-") == TRUE),TRUE)


  expect_equal(nrow(test.nosoiA$table.hosts_A), nrow(test.nosoiA$table.state_A))

  expect_equal(test.nosoiA$total.time, 26)

  expect_equal(test.nosoiA$N.infected_A, 627)
  expect_equal(test.nosoiA$N.infected_B, 520)

  expect_equal(test.nosoiA$type, "dualContinuous")

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
  test.nosoiA <- nosoiSim(type="dual",structure=TRUE, continuous = TRUE,
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
                          prefix.host.B="V",

                          progress.bar=TRUE,
                          print.step=10)

  full.results.nosoi <- rbindlist(list(test.nosoiA$table.hosts_A,test.nosoiA$table.hosts_B))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$table.state_A,test.nosoiA$table.state_B))

  expect_equal(nrow(test.nosoiA$table.hosts_A), nrow(test.nosoiA$table.state_A))

  expect_equal(test.nosoiA$total.time, 4)

  expect_equal(test.nosoiA$N.infected_A, 0)

  expect_equal(test.nosoiA$N.infected_B, 1)

  expect_equal(test.nosoiA$type, "dualContinuous")

  #Movement

  H1_moves <- subset(full.results.nosoi.state, hosts.ID == "V-1")
  expect_equal(nrow(H1_moves),2)
})
