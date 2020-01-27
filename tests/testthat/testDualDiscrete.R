context("Testing dual-host with discrete structure")

test_that("Movement is coherent with single introduction, constant pMove, diff pExit", {
  skip_if_not_installed("igraph")
  library(igraph)

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

  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts,test.nosoiA$host.info.B$table.hosts))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 7)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 20)

  expect_equal(test.nosoiA$host.info.A$N.infected, 75)
  expect_equal(test.nosoiA$host.info.B$N.infected, 67)

  expect_equal(test.nosoiA$type, "dual")
  expect_equal(test.nosoiA$host.info.A$popStructure, "discrete")
  expect_equal(test.nosoiA$host.info.B$popStructure, "discrete")

  #Movement
  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "H-3")),2)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "H-3")$state,c("A","C"))

  Where.when.exit = subset(full.results.nosoi,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.exit, current.in == "A")$N,integer(0))
  expect_equal(subset(Where.when.exit, current.in == "B")$N,21)
  expect_equal(subset(Where.when.exit, current.in == "C")$N,29)

  #Test output

  test <- summary(test.nosoiA)

  expect_equal(test$R0$N.inactive.A, 26)
  expect_equal(test$dynamics[21]$t, 11)
  expect_equal(test$dynamics[21]$Count, 2)
  expect_equal(test$dynamics[21]$type, "V")
  expect_equal(test$dynamics[21]$state, "A")

  #Get host table
  test.hostTable.A <- getTableHosts(test.nosoiA, pop="A")
  expect_equal(test.hostTable.A[35]$inf.by, "V-11")

  test.hostTable.B <- getTableHosts(test.nosoiA, pop="B")
  expect_equal(test.hostTable.B[35]$inf.by, "H-16")

  #Get state table
  test.stateTable.A <- getTableState(test.nosoiA, pop="A")
  expect_equal(test.stateTable.A[52]$hosts.ID, "H-40")

  test.stateTable.B <- getTableState(test.nosoiA, pop="B")
  expect_equal(test.stateTable.B[52]$hosts.ID, "V-44")

  expect_error(test.stateTable.A <- getTableHosts(test.nosoiA, pop="Z"),
               "Population Z is not recognized.")

})

test_that("Transmission is coherent with single introduction (host A) differential according to host, shared parameter", {
  skip_if_not_installed("igraph")
  library(igraph)

  #Host A

  t_infectA_fct <- function(x){rnorm(x,mean = 12,sd=3)}
  pTrans_hostA <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=sin(pi*t/t_infectA)}
    if(t/t_infectA > 1){p=0}
    return(p)
  }

  p_Move_fctA  <- function(t){return(0.1)}

  p_Exit_fctA  <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=0}
    if(t/t_infectA > 1){p=1}
    return(p)
  }

  time_contact_A = function(t){sample(c(0,1,2),1,prob=c(0.2,0.4,0.4))}

  #Host B
  t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_Exit_fct_B  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0.2)}
    if(current.in=="C"){return(1)}}

  pTrans_hostB <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(6262)
  test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                          length.sim=40,
                          max.infected.A=100,
                          max.infected.B=200,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A="A",
                          init.structure.B=NA,
                          structure.matrix.A=transition.matrix,
                          structure.matrix.B=transition.matrix,

                          pExit.A = p_Exit_fctA,
                          param.pExit.A = list(t_infectA = t_infectA_fct),
                          pMove.A=p_Move_fctA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          timeDep.pExit.A=FALSE,
                          nContact.A = time_contact_A,
                          param.nContact.A = NA,
                          timeDep.nContact.A=FALSE,
                          pTrans.A = pTrans_hostA,
                          param.pTrans.A = list(t_infectA=t_infectA_fct),
                          timeDep.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B = p_Exit_fct_B,
                          param.pExit.B = NA,
                          timeDep.pExit.B=FALSE,
                          diff.pExit.B=TRUE,
                          pMove.B=NA,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          nContact.B = time_contact_B,
                          param.nContact.B = NA,
                          timeDep.nContact.B=FALSE,
                          pTrans.B = pTrans_hostB,
                          param.pTrans.B = list(p_max=p_max_fct_B,
                                                t_incub=t_incub_fct_B),
                          timeDep.pTrans.B=FALSE,
                          prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts[,c(1:7)],test.nosoiA$host.info.B$table.hosts[,c(1:7)]))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 9)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 27)

  expect_equal(test.nosoiA$host.info.A$N.infected, 109)
  expect_equal(test.nosoiA$host.info.B$N.infected, 223)

  expect_equal(colnames(test.nosoiA$host.info.A$table.hosts)[10],"t_infectA")

  #Movement
  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "H-3")),2)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "H-3")$state,c("B","A"))

  Where.when.exit.B <- subset(test.nosoiA$host.info.B$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  Where.when.exit.A <- subset(test.nosoiA$host.info.A$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.exit.A, current.in == "A")$N,2)
  expect_equal(subset(Where.when.exit.B, current.in == "B")$N,40)

  test.B <- test.nosoiA$host.info.B$table.state %>% group_by(hosts.ID) %>% summarise(N=length(hosts.ID))
  expect_equal(unique(test.B$N), 1) #B should not be moving

  test.A <- test.nosoiA$host.info.A$table.state %>% group_by(hosts.ID) %>% summarise(N=length(hosts.ID))
  expect_equal(length(test.A$N) > 1, TRUE)
})

test_that("Transmission is coherent with single introduction (host A) differential according to host, shared parameter, A does not move", {
  skip_if_not_installed("igraph")
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

  p_Exit_fct_B  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0.2)}
    if(current.in=="C"){return(1)}}

  p_Move_fctB  <- function(t){return(0.1)}

  pTrans_hostB <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(6262)
  test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                          length.sim=40,
                          max.infected.A=100,
                          max.infected.B=200,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A="A",
                          init.structure.B=NA,
                          structure.matrix.A=transition.matrix,
                          structure.matrix.B=transition.matrix,

                          pExit.A = p_Exit_fctA,
                          param.pExit.A = list(t_infectA = t_infectA_fct),
                          pMove.A=NA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          timeDep.pExit.A=FALSE,
                          nContact.A = time_contact_A,
                          param.nContact.A = NA,
                          timeDep.nContact.A=FALSE,
                          pTrans.A = pTrans_hostA,
                          param.pTrans.A = list(t_infectA=t_infectA_fct),
                          timeDep.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B = p_Exit_fct_B,
                          param.pExit.B = NA,
                          timeDep.pExit.B=FALSE,
                          diff.pExit.B=TRUE,
                          pMove.B=p_Move_fctB,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          nContact.B = time_contact_B,
                          param.nContact.B = NA,
                          timeDep.nContact.B=FALSE,
                          pTrans.B = pTrans_hostB,
                          param.pTrans.B = list(p_max=p_max_fct_B,
                                                t_incub=t_incub_fct_B),
                          timeDep.pTrans.B=FALSE,
                          prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts[,c(1:7)],test.nosoiA$host.info.B$table.hosts[,c(1:7)]))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 11)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 21)

  expect_equal(test.nosoiA$host.info.A$N.infected, 125)
  expect_equal(test.nosoiA$host.info.B$N.infected, 238)

  expect_equal(colnames(test.nosoiA$host.info.A$table.hosts)[10],"t_infectA")

  #Movement
  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "H-3")),1)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "H-3")$state,c("A"))

  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "V-2")),2)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "V-2")$state,c("A","C"))

  Where.when.exit.B <- subset(test.nosoiA$host.info.B$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  Where.when.exit.A <- subset(test.nosoiA$host.info.A$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.exit.A, current.in == "A")$N,3)
  expect_equal(subset(Where.when.exit.B, current.in == "B")$N,25)

  expect_error(
    test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                            length.sim=40,
                            max.infected.A=100,
                            max.infected.B=200,
                            init.individuals.A=1,
                            init.individuals.B=0,
                            init.structure.A="A",
                            init.structure.B=NA,
                            structure.matrix.A=transition.matrix,
                            structure.matrix.B=transition.matrix,

                            pExit.A = p_Exit_fctA,
                            param.pExit.A = list(t_infectA = t_infectA_fct),
                            pMove.A=NA,
                            param.pMove.A=NA,
                            timeDep.pMove.A=FALSE,
                            diff.pMove.A=FALSE,
                            timeDep.pExit.A=FALSE,
                            nContact.A = time_contact_A,
                            param.nContact.A = NA,
                            timeDep.nContact.A=FALSE,
                            pTrans.A = pTrans_hostA,
                            param.pTrans.A = list(t_infectA=t_infectA_fct),
                            timeDep.pTrans.A=FALSE,
                            prefix.host.A="H",

                            pExit.B = p_Exit_fct_B,
                            param.pExit.B = NA,
                            timeDep.pExit.B=FALSE,
                            diff.pExit.B=TRUE,
                            pMove.B=NA,
                            param.pMove.B=NA,
                            timeDep.pMove.B=FALSE,
                            diff.pMove.B=FALSE,
                            nContact.B = time_contact_B,
                            param.nContact.B = NA,
                            timeDep.nContact.B=FALSE,
                            pTrans.B = pTrans_hostB,
                            param.pTrans.B = list(p_max=p_max_fct_B,
                                                  t_incub=t_incub_fct_B),
                            timeDep.pTrans.B=FALSE,
                            prefix.host.B="V"),
    "At least one host must move."
  )
})

test_that("Epidemics dying out", {
  skip_if_not_installed("igraph")
  library(igraph)

  #Host A

  t_infectA_fct <- function(x){rnorm(x,mean = 12,sd=3)}
  pTrans_hostA <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=sin(pi*t/t_infectA)}
    if(t/t_infectA > 1){p=0}
    return(p)
  }

  p_Move_fctA  <- function(t){return(0.1)}

  p_Exit_fctA  <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=0}
    if(t/t_infectA > 1){p=1}
    return(p)
  }

  time_contact_A = function(t){sample(c(0,1,2),1,prob=c(0.2,0.4,0.4))}

  #Host B
  t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_Exit_fct_B  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0.2)}
    if(current.in=="C"){return(1)}}

  pTrans_hostB <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(101)
  test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                          length.sim=40,
                          max.infected.A=100,
                          max.infected.B=200,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A="A",
                          init.structure.B=NA,
                          structure.matrix.A=transition.matrix,
                          structure.matrix.B=transition.matrix,

                          pExit.A = p_Exit_fctA,
                          param.pExit.A = list(t_infectA = t_infectA_fct),
                          pMove.A=p_Move_fctA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          timeDep.pExit.A=FALSE,
                          nContact.A = time_contact_A,
                          param.nContact.A = NA,
                          timeDep.nContact.A=FALSE,
                          pTrans.A = pTrans_hostA,
                          param.pTrans.A = list(t_infectA=t_infectA_fct),
                          timeDep.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B = p_Exit_fct_B,
                          param.pExit.B = NA,
                          timeDep.pExit.B=FALSE,
                          diff.pExit.B=TRUE,
                          pMove.B=NA,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          nContact.B = time_contact_B,
                          param.nContact.B = NA,
                          timeDep.nContact.B=FALSE,
                          pTrans.B = pTrans_hostB,
                          param.pTrans.B = list(p_max=p_max_fct_B,
                                                t_incub=t_incub_fct_B),
                          timeDep.pTrans.B=FALSE,
                          prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts[,c(1:7)],test.nosoiA$host.info.B$table.hosts[,c(1:7)]))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 2)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 12)

  expect_equal(test.nosoiA$host.info.A$N.infected, 1)
  expect_equal(test.nosoiA$host.info.B$N.infected, 8)

  #Movement
  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "H-1")),3)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "H-1")$state,c("A","B","C"))
})

test_that("start with host B", {
  skip_if_not_installed("igraph")
  library(igraph)

  #Host B

  t_infectA_fct <- function(x){rnorm(x,mean = 12,sd=3)}
  pTrans_hostA <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=sin(pi*t/t_infectA)}
    if(t/t_infectA > 1){p=0}
    return(p)
  }

  p_Move_fctA  <- function(t){return(0.1)}

  p_Exit_fctA  <- function(t,t_infectA){
    if(t/t_infectA <= 1){p=0}
    if(t/t_infectA > 1){p=1}
    return(p)
  }

  time_contact_A = function(t){sample(c(0,1,2),1,prob=c(0.2,0.4,0.4))}

  #Host B
  t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_Exit_fct_B  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0.2)}
    if(current.in=="C"){return(1)}}

  pTrans_hostB <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact_B = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(101)
  test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                          length.sim=40,
                          max.infected.A=100,
                          max.infected.B=200,
                          init.individuals.A=0,
                          init.individuals.B=1,
                          init.structure.A=NA,
                          init.structure.B="A",
                          structure.matrix.A=transition.matrix,
                          structure.matrix.B=transition.matrix,

                          pExit.A = p_Exit_fctA,
                          param.pExit.A = list(t_infectA = t_infectA_fct),
                          pMove.A=p_Move_fctA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          timeDep.pExit.A=FALSE,
                          nContact.A = time_contact_A,
                          param.nContact.A = NA,
                          timeDep.nContact.A=FALSE,
                          pTrans.A = pTrans_hostA,
                          param.pTrans.A = list(t_infectA=t_infectA_fct),
                          timeDep.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B = p_Exit_fct_B,
                          param.pExit.B = NA,
                          timeDep.pExit.B=FALSE,
                          diff.pExit.B=TRUE,
                          pMove.B=NA,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          nContact.B = time_contact_B,
                          param.nContact.B = NA,
                          timeDep.nContact.B=FALSE,
                          pTrans.B = pTrans_hostB,
                          param.pTrans.B = list(p_max=p_max_fct_B,
                                                t_incub=t_incub_fct_B),
                          timeDep.pTrans.B=FALSE,
                          prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts[,c(1:7)],test.nosoiA$host.info.B$table.hosts[,c(1:7)]))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 9)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 17)

  expect_equal(test.nosoiA$host.info.A$N.infected, 118)
  expect_equal(test.nosoiA$host.info.B$N.infected, 216)

  #Movement
  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "H-1")),3)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "H-1")$state,c("A","B","C"))
})

test_that("Transmission is coherent with single introduction (host A) differential according to host, shared parameter, A does not move", {
  skip_if_not_installed("igraph")
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

  time_contact_A = function(t){round(rnorm(1, 3, 1), 0)}

  #Host B
  t_incub_fct_B <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct_B <- function(x){rbeta(x,shape1 = 5,shape2=2)}

  p_Exit_fct_B  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0.2)}
    if(current.in=="C"){return(1)}}

  p_Move_fctB  <- function(t){return(0.1)}

  pTrans_hostB <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact_B <- function(t, current.in, host.count.A){

    temp.val = 30 - host.count.A

    if(temp.val <= 0) {
      return(0)
    }
    if(temp.val >= 0) {
      if(current.in=="A"){
        return(round((temp.val/30)*rnorm(1, 3, 1), 0))}
      if(current.in=="B"){return(0)}
      if(current.in=="C"){
        return(round((temp.val/30)*rnorm(1,1, 1), 0))}
    }
  }

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(6262)
  test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                          length.sim=20,
                          max.infected.A=100,
                          max.infected.B=200,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A="A",
                          init.structure.B=NA,
                          structure.matrix.A=transition.matrix,
                          structure.matrix.B=transition.matrix,

                          pExit.A = p_Exit_fctA,
                          param.pExit.A = list(t_infectA = t_infectA_fct),
                          pMove.A=NA,
                          param.pMove.A=NA,
                          timeDep.pMove.A=FALSE,
                          diff.pMove.A=FALSE,
                          timeDep.pExit.A=FALSE,

                          nContact.A = time_contact_A,
                          param.nContact.A = NA,

                          pTrans.A = pTrans_hostA,
                          param.pTrans.A = list(t_infectA=t_infectA_fct),
                          timeDep.pTrans.A=FALSE,
                          prefix.host.A="H",

                          pExit.B = p_Exit_fct_B,
                          param.pExit.B = NA,
                          timeDep.pExit.B=FALSE,
                          diff.pExit.B=TRUE,
                          pMove.B=p_Move_fctB,
                          param.pMove.B=NA,
                          timeDep.pMove.B=FALSE,
                          diff.pMove.B=FALSE,
                          nContact.B = time_contact_B,
                          param.nContact.B = NA,
                          timeDep.nContact.B=FALSE,
                          diff.nContact.B=TRUE,
                          hostCount.nContact.B=TRUE,
                          pTrans.B = pTrans_hostB,
                          param.pTrans.B = list(p_max=p_max_fct_B,
                                                t_incub=t_incub_fct_B),
                          timeDep.pTrans.B=FALSE,
                          prefix.host.B="V")


  full.results.nosoi <- rbindlist(list(test.nosoiA$host.info.A$table.hosts[,c(1:7)],test.nosoiA$host.info.B$table.hosts[,c(1:7)]))
  full.results.nosoi.state <- rbindlist(list(test.nosoiA$host.info.A$table.state,test.nosoiA$host.info.B$table.state))

  g <- graph.data.frame(full.results.nosoi[inf.by != "NA-1",c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 7)

  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts$inf.by,"H-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.A$table.hosts[-1]$inf.by,"V-") == TRUE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts$inf.by,"V-") == FALSE),TRUE)
  expect_equal(all(str_detect(test.nosoiA$host.info.B$table.hosts[-1]$inf.by,"H-") == TRUE),TRUE)

  expect_equal(test.nosoiA$total.time, 14)

  expect_equal(test.nosoiA$host.info.A$N.infected, 31)
  expect_equal(test.nosoiA$host.info.B$N.infected, 217)

  expect_equal(colnames(test.nosoiA$host.info.A$table.hosts)[10],"t_infectA")

  #Movement
  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "H-3")),1)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "H-3")$state,c("A"))

  expect_equal(nrow(subset(full.results.nosoi.state, hosts.ID == "V-2")),2)
  expect_equal(subset(full.results.nosoi.state, hosts.ID == "V-2")$state,c("A","B"))

  Where.when.exit.B <- subset(test.nosoiA$host.info.B$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  Where.when.exit.A <- subset(test.nosoiA$host.info.A$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.exit.A, current.in == "A")$N,1)
  expect_equal(subset(Where.when.exit.B, current.in == "B")$N,5)

  #Number of host at each loc at each time
  out = getDynamic(test.nosoiA) #  ggplot(out, aes(x=t,y=Count,lty=state,color=type)) + geom_line() + geom_hline(yintercept=30)
  expect_equal(subset(out,state=="A" & t==15 & type=="H")$Count,28)
  expect_equal(length(subset(out,state=="B" & t==15 & type=="H")$Count),0)
  expect_equal(subset(out,state=="C" & t==15 & type=="H")$Count,2)

})
