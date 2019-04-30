context("Testing dual-host without structure")

test_that("Transmission is coherent with single introduction (host A) same for both hosts", {
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

  set.seed(66)
  test.nosoiA <- nosoiSim(type="dual",structure=FALSE,
                       length.sim=40,
                       max.infected.A=100,
                       max.infected.B=100,
                       init.individuals.A=1,
                       init.individuals.B=0,

                       pExit.A = p_Exit_fct,
                       param.pExit.A = NA,
                       timeDep.pExit.A=FALSE,
                       timeContact.A = time_contact,
                       param.timeContact.A = NA,
                       timeDep.timeContact.A=FALSE,
                       pTrans.A = proba,
                       param.pTrans.A = list(p_max=p_max_fct,
                                             t_incub=t_incub_fct),
                       timeDep.pTrans.A=FALSE,
                       prefix.host.A="H",

                       pExit.B = p_Exit_fct,
                       param.pExit.B = NA,
                       timeDep.pExit.B=FALSE,
                       timeContact.B = time_contact,
                       param.timeContact.B = NA,
                       timeDep.timeContact.B=FALSE,
                       pTrans.B = proba,
                       param.pTrans.B = list(p_max=p_max_fct,
                                             t_incub=t_incub_fct),
                       timeDep.pTrans.B=FALSE,
                       prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$table.hosts_A,test.nosoiA$table.hosts_B))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

  expect_equal(all(str_detect(test.nosoiA$table.hosts_A$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_A[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 20)

  expect_equal(test.nosoiA$N.infected_A, 126)
  expect_equal(test.nosoiA$N.infected_B, 87)

  expect_equal(test.nosoiA$type, "dualNone")
})

test_that("Transmission is coherent with single introduction (host A) differential according to host, shared parameter", {
  library(igraph)

  #Host A

  t_infectA_fct <- function(x){rnorm(x,mean = 12,sd=3)}
  pTrans_hostA <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=sin(pi*t/t_infectA)}
    if(t/t_infectA > 1){p=0}
    return(p)
  }

  p_Exit_fctA  <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=0}
    if(t/t_infectA > 1){p=1}
    return(p)
  }

  time_contact_A = function(t){sample(c(0,1,2),1,prob=c(0.2,0.4,0.4))}

  #Host B
  t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_Exit_fct_B  <- function(t){return(0.08)}

  pTrans_hostB <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(150)
  test.nosoiA <- nosoiSim(type="dual",structure=FALSE,
                          length.sim=40,
                          max.infected.A=100,
                          max.infected.B=200,
                          init.individuals.A=1,
                          init.individuals.B=0,

                          pExit.A = p_Exit_fctA,
                          param.pExit.A = list(t_infectA = t_infectA_fct),
                          timeDep.pExit.A=FALSE,
                          timeContact.A = time_contact_A,
                          param.timeContact.A = NA,
                          timeDep.timeContact.A=FALSE,
                          pTrans.A = pTrans_hostA,
                          param.pTrans.A = list(t_infectA=t_infectA_fct),
                          timeDep.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B = p_Exit_fct_B,
                          param.pExit.B = NA,
                          timeDep.pExit.B=FALSE,
                          timeContact.B = time_contact_B,
                          param.timeContact.B = NA,
                          timeDep.timeContact.B=FALSE,
                          pTrans.B = pTrans_hostB,
                          param.pTrans.B = list(p_max=p_max_fct_B,
                                                t_incub=t_incub_fct_B),
                          timeDep.pTrans.B=FALSE,
                          prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$table.hosts_A[,c(1,2)],test.nosoiA$table.hosts_B[,c(1,2)]))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 10)

  expect_equal(all(str_detect(test.nosoiA$table.hosts_A$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_A[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 17)

  expect_equal(test.nosoiA$N.infected_A, 105)
  expect_equal(test.nosoiA$N.infected_B, 226)

  expect_equal(colnames(test.nosoiA$table.hosts_A)[6],"t_infectA")
})

test_that("Transmission is coherent with single introduction (host A) differential according to host, shared parameter, time dependancy for host B pExit", {
  library(igraph)

  #Host A

  t_infectA_fct <- function(x){rnorm(x,mean = 12,sd=3)}
  pTrans_hostA <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=sin(pi*t/t_infectA)}
    if(t/t_infectA > 1){p=0}
    return(p)
  }

  p_Exit_fctA  <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=0}
    if(t/t_infectA > 1){p=1}
    return(p)
  }

  time_contact_A = function(t){sample(c(0,1,2),1,prob=c(0.2,0.4,0.4))}

  #Host B
  t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_Exit_fct_B  <- function(t,prestime){(sin(prestime/12)+1)/5}

  pTrans_hostB <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}

  set.seed(90)
  test.nosoiA <- nosoiSim(type="dual",structure=FALSE,
                          length.sim=40,
                          max.infected.A=100,
                          max.infected.B=200,
                          init.individuals.A=1,
                          init.individuals.B=0,

                          pExit.A = p_Exit_fctA,
                          param.pExit.A = list(t_infectA = t_infectA_fct),
                          timeDep.pExit.A=FALSE,
                          timeContact.A = time_contact_A,
                          param.timeContact.A = NA,
                          timeDep.timeContact.A=FALSE,
                          pTrans.A = pTrans_hostA,
                          param.pTrans.A = list(t_infectA=t_infectA_fct),
                          timeDep.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B = p_Exit_fct_B,
                          param.pExit.B = NA,
                          timeDep.pExit.B=TRUE,
                          timeContact.B = time_contact_B,
                          param.timeContact.B = NA,
                          timeDep.timeContact.B=FALSE,
                          pTrans.B = pTrans_hostB,
                          param.pTrans.B = list(p_max=p_max_fct_B,
                                                t_incub=t_incub_fct_B),
                          timeDep.pTrans.B=FALSE,
                          prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$table.hosts_A[,c(1,2)],test.nosoiA$table.hosts_B[,c(1,2)]))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 12)

  expect_equal(all(str_detect(test.nosoiA$table.hosts_A$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_A[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$table.hosts_B[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 39)

  expect_equal(test.nosoiA$N.infected_A, 71)
  expect_equal(test.nosoiA$N.infected_B, 221)

})