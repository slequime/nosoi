context("Testing single-host without structure")

test_that("Transmission is coherent with single introduction, constant pExit and pTrans", {
  skip_if_not_installed("igraph")
  library(igraph)

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t){return(0.2)}

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(805)
  test.nosoiA <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=1,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = NA,
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )
  g <- graph.data.frame(getHostData(test.nosoiA, "table.host")[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 12)

})


test_that("Transmission is coherent with single introduction, constant pExit and pTrans", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(805)
  test.nosoiA <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=1,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )
  g <- graph.data.frame(getHostData(test.nosoiA, "table.host")[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

  test <- summary(test.nosoiA)

  expect_equal(test$R0$N.inactive, 16)
  expect_equal(test$R0$R0.mean, 0)
  expect_equal(test$dynamics[12]$t, 11)
  expect_equal(test$dynamics[12]$Count, 39)
  expect_equal(test$dynamics[12]$type, "H")
  expect_equal(test$cumulative[14]$t, 13)
  expect_equal(test$cumulative[14]$Count, 89)
  expect_equal(test$cumulative[14]$type, "H")

  #get Tables
  #Get host table
  test.hostTable.A <- getTableHosts(test.nosoiA)
  expect_equal(test.hostTable.A[35]$out.time, 10)

  #Check errors
  expect_error(test.stateTable.A <- getTableState(test.nosoiA),
               "There is no state information kept when the host population A has no structure.")

  expect_error(test.stateTable.A <- getTableHosts(test.nosoiA, pop="B"),
               "There are no other hosts than 'A' in a single-host simulation.")

})

test_that("Transmission is coherent with single introduction, simple pExit and pTrans", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){plogis(t,20,2)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(805)
  test.nosoiB <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=1,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )
  g <- graph.data.frame(getHostData(test.nosoiB, "table.host")[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)
})

test_that("Transmission is coherent with single introduction, complex pExit and pTrans", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_exit_param1 <- function(x){rnorm(x,mean = 10,sd=2)}

  p_Exit_fct  <- function(t,pExit.param1){plogis(t,pExit.param1,2)}
  # p_exit_bis <- function(x){p_Exit_fct(t=x[1], pExit.param1=x[2])}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(805)
  test.nosoiC <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=1,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit = list(pExit.param1=p_exit_param1)
  )
  g <- graph.data.frame(getHostData(test.nosoiC, "table.host")[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)
})

test_that("Transmission is coherent with multiple introductions, constant pExit and pTrans", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  # p_Exit_fct  <- function(x){rep(0.08,length(x))}
  p_Exit_fct  <- function(t){return(0.08)}


  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(805)
  test.nosoiA <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=3,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  g <- graph.data.frame(getHostData(test.nosoiA, "table.host")[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 3)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

})

test_that("Transmission is coherent with multiple introductions, simple pExit and pTrans", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){plogis(t,20,2)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(805)
  test.nosoiB <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=3,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )
  g <- graph.data.frame(getHostData(test.nosoiB, "table.host")[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 3)
  expect_equal(diameter(g, directed=F, weights=NA), 4)
})

test_that("Transmission is coherent with multiple introductions, complex pExit and pTrans", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_exit_param1 <- function(x){rnorm(x,mean = 10,sd=2)}

  p_Exit_fct  <- function(t,pExit.param1){plogis(t,pExit.param1,2)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(805)
  test.nosoiC <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=3,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit = list(pExit.param1=p_exit_param1)
  )
  g <- graph.data.frame(getHostData(test.nosoiC, "table.host")[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 3)
  expect_equal(diameter(g, directed=F, weights=NA), 6)
})

test_that("Dying out epidemic", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(111)
  test.nosoiA <- nosoiSim(type="single", popStructure="none",
                          length=40,
                          max.infected=100,
                          init.individuals=1,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  expect_equal(nrow(getHostData(test.nosoiA, "table.host")),1)
  expect_equal(test.nosoiA$total.time,4)
})
