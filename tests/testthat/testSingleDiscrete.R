context("Testing single-host with discrete structure")

test_that("Error message pops up when structure.matrix or init.structure is badly formed", {
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

  # transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))
  transition.matrix = c("A","B")

  expect_error(
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
    ),
    "structure.matrix should be a matrix."
  )

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0,0.8,0.1,0.3),nrow = 3, ncol = 4,dimnames=list(c("A","B","C"),c("A","B","C","D")))

  expect_error(
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
    ),
    "structure.matrix should have the same number of rows and columns."
  )

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","D")))

  expect_error(
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
    ),
    "structure.matrix rows and columns should have the same names."
  )

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.8,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  expect_error(
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
    ),
    "structure.matrix rows should sum up to 1."
  )

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                            length=20,
                            max.infected=100,
                            init.individuals=1,
                            init.structure="D",
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
    ),
    "init.structure should be a state present in structure.matrix."
  )

})

test_that("Movement is coherent with single introduction, constant pMove", {
  skip_if_not_installed("igraph")
  library(igraph)

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

  ## Output
  expect_output(print(test.nosoiA), "a single host with a discrete structure")

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")),3)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")$state,c("A","C","B"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.at.end = getHostData(test.nosoiA, "table.hosts") %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.at.end, current.in == "A")$N,62)
  expect_equal(subset(Where.at.end, current.in == "B")$N,17)
  expect_equal(subset(Where.at.end, current.in == "C")$N,43)

  #Summary
  test <- summary(test.nosoiA)

  expect_equal(test$R0$N.inactive, 31)
  expect_equal(test$dynamics[21]$t, 14)
  expect_equal(test$dynamics[21]$Count, 11)
  expect_equal(test$dynamics[21]$type, "H")
  expect_equal(test$dynamics[21]$state, "C")

  #Get host table
  test.hostTable.A <- getTableHosts(test.nosoiA)
  expect_equal(test.hostTable.A[35]$out.time, 14)

  #Get state table
  test.stateTable.A <- getTableState(test.nosoiA)
  expect_equal(test.stateTable.A[52]$hosts.ID, "H-26")

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})

#--------------------------------------------------------------------------------------------------------------------------------------------

test_that("Movement is coherent with single introduction, complex pMove", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}
  p_Move_param1_fct <- function(x){rnorm(x,mean = 10,sd=2)}

  p_Move_fct  <- function(t,pMove.param1){plogis(t,pMove.param1,2)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(750)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=20,
                          max.infected=1000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          pMove=p_Move_fct,
                          param.pMove=list(pMove.param1=p_Move_param1_fct),
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),11)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")$state,c("A","C","A","C","A","B","C","B","C","A","C"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.at.end = getHostData(test.nosoiA, "table.hosts") %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.at.end, current.in == "A")$N,83)
  expect_equal(subset(Where.at.end, current.in == "B")$N,19)
  expect_equal(subset(Where.at.end, current.in == "C")$N,16)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})

test_that("Movement is coherent with single introduction, constant but different pMove, 1 loc (C) is sink. Ce tombeau sera votre tombeau !", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}

  p_Move_fct  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0.3)}
    if(current.in=="C"){return(0)}}
  #

  # p_Move_param1_fct <- function(x){rnorm(x,mean = 10,sd=2)}
  #
  # p_Move_fct  <- function(t,current.in,pMove.param1){
  #   if(current.in=="A"){return(plogis(t,1+pMove.param1,2))}
  #   if(current.in=="B"){return(plogis(t,2+pMove.param1,2))}
  #   if(current.in=="C"){return(plogis(t,pMove.param1,2)/1000)}}

  # param.pMove = list(pMove.param1=p_Move_param1_fct)
  param.pMove = NA

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(750)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=20,
                          max.infected=10000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          diff.pMove=TRUE,
                          pMove=p_Move_fct,
                          param.pMove=param.pMove,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),2)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")$state,c("A","C"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.at.end = getHostData(test.nosoiA, "table.hosts") %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  expect_equal(subset(Where.at.end, current.in == "C")$N,138)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))

})

test_that("Error message pops up if different pMove poorly formatted", {
  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}


  p_Move_fct  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="C"){return(0)}}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                            length=20,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            diff.pMove=TRUE,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            pExit=p_Exit_fct,
                            param.pExit=NA
    ),
    "pMove should have a realisation for each possible state. diff.pMove is TRUE."
  )
})

test_that("Movement is coherent with single introduction, complex and different pMove, 1 loc (C) is sink.", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}

  p_Move_param1_fct <- function(x){rnorm(x,mean = 10,sd=2)}
  p_Move_fct  <- function(t,current.in,pMove.param1){
    if(current.in=="A"){return(plogis(t,1+pMove.param1,2))}
    if(current.in=="B"){return(plogis(t,2+pMove.param1,2))}
    if(current.in=="C"){return(plogis(t,pMove.param1,2)/1000)}}
  param.pMove = list(pMove.param1=p_Move_param1_fct)

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(750)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=20,
                          max.infected=10000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          diff.pMove=TRUE,
                          pMove=p_Move_fct,
                          param.pMove=param.pMove,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),2)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")$state,c("A","C"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.at.end = getHostData(test.nosoiA, "table.hosts") %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  expect_equal(subset(Where.at.end, current.in == "C")$N,52)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})

test_that("Error message pops out when missing state in diff functions", {

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  p_Exit_fct  <- function(t,current.in){
    if(current.in=="A"){return(0)}
    if(current.in=="C"){return(1)}}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                            length=20,
                            max.infected=1000,
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
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA
    ),
    "pExit should have a realisation for each possible state. diff.pExit is TRUE."
  )

  p_Exit_fct  <- function(x){return(0.08)}

  time_contact <- function(t,current.in){
    if(current.in=="A"){
      return(round(rnorm(1, 3, 1), 0))}
    if(current.in=="C"){
      return(round(rnorm(1, 6, 1), 0))}
  }

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                            length=20,
                            max.infected=1000,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            diff.nContact=TRUE,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=NA,
                            pExit=p_Exit_fct,
                            param.pExit=NA
    ),
    "nContact should have a realisation for each possible state. diff.nContact is TRUE."
  )

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  proba <- function(t,current.in,p_max,t_incub){
    if(current.in=="A"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=p_max}
      return(p)}
    if(current.in=="C"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=p_max}
      return(p)}
  }

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                            length=20,
                            max.infected=1000,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            nContact=time_contact,
                            param.nContact=NA,
                            diff.pTrans=TRUE,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            pExit=p_Exit_fct,
                            param.pExit=NA
    ),
    "pTrans should have a realisation for each possible state. diff.pTrans is TRUE."
  )

})


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
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=20,
                          max.infected=1000,
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
                          diff.pExit=TRUE,
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)], directed = F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 7)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")),2)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")$state,c("A","C"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.when.exit = subset(getHostData(test.nosoiA, "table.hosts"),active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.exit, current.in == "A")$N,integer(0))
  expect_equal(subset(Where.when.exit, current.in == "B")$N,44)
  expect_equal(subset(Where.when.exit, current.in == "C")$N,34)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})

test_that("Movement is coherent with single introduction, constant pMove, diff pTrans ", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,current.in,p_max,t_incub){
    if(current.in=="A"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=p_max}
      return(p)}
    if(current.in=="B"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=0}
      return(p)}
    if(current.in=="C"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=p_max}
      return(p)}
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(805)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=20,
                          max.infected=1000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          nContact=time_contact,
                          param.nContact=NA,
                          diff.pTrans = TRUE,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 7)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")),2)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")$state,c("A","C"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.when.infected = subset(getHostData(test.nosoiA, "table.hosts"),active==0) %>% group_by(inf.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.infected, inf.in == "B")$N,integer(0))
  expect_equal(subset(Where.when.infected, inf.in == "A")$N,56)
  expect_equal(subset(Where.when.infected, inf.in == "C")$N,17)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))

})

test_that("Movement is coherent with single introduction, constant pMove, diff nContact ", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact <- function(t,current.in){
    if(current.in=="A"){
      return(round(rnorm(1, 3, 1), 0))}
    if(current.in=="B"){return(0)}
    if(current.in=="C"){
      return(round(rnorm(1, 6, 1), 0))}
  }

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(805)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=20,
                          max.infected=1000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.nContact=TRUE,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 8)

  expect_equal(nrow(getHostData(test.nosoiA, "table.hosts")),612)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")),2)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-3")$state,c("A","C"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.when.infected = subset(getHostData(test.nosoiA, "table.hosts"),active==0) %>% group_by(inf.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.infected, inf.in == "B")$N,integer(0))
  expect_equal(subset(Where.when.infected, inf.in == "A")$N,27)
  expect_equal(subset(Where.when.infected, inf.in == "C")$N,67)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})

test_that("Movement is coherent with single introduction, all parameters are diff", {
  skip_if_not_installed("igraph")
  library(igraph)

  p_Move_fct  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0)}
    if(current.in=="C"){return(0.2)}
    if(current.in=="D"){return(0.3)}
  }

  p_Exit_fct  <- function(t,current.in){
    if(current.in=="A"){return(0.1)}
    if(current.in=="B"){return(0.1)}
    if(current.in=="C"){return(0)}
    if(current.in=="D"){return(1)}
  }

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  proba <- function(t,current.in,p_max,t_incub){
    if(current.in=="A"){
      if(t <= t_incub-3){p=0}
      if(t >= t_incub){p=p_max}
      return(p)}
    if(current.in=="B"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=0}
      return(p)}
    if(current.in=="C"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=1}
      return(p)}
    if(current.in=="D"){
      if(t <= t_incub){p=0}
      if(t >= t_incub){p=p_max}
      return(p)}
  }

  time_contact <- function(t,current.in){
    if(current.in=="A"){
      return(0)}
    if(current.in=="B"){return(0)}
    if(current.in=="C"){
      return(round(rnorm(1, 2, 1), 0))}
    if(current.in=="D"){
      return(round(rnorm(1, 3, 1), 0))}
  }

  transition.matrix = matrix(c(0,0.1,0.6,0.3,0.1,0,0.1,0.3,0.6,0.5,0,0.4,0.3,0.4,0.3,0),nrow = 4, ncol = 4,dimnames=list(c("A","B","C","D"),c("A","B","C","D")))
  # transition.matrix = matrix(c(0,2,3,4,5,0,7,8,9,10,0,12,13,14,15,0),nrow = 4, ncol = 4,dimnames=list(c("A","B","C","D"),c("A","B","C","D")))

  set.seed(186)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=40,
                          max.infected=1000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          diff.pMove=TRUE,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.nContact=TRUE,
                          nContact=time_contact,
                          param.nContact=NA,
                          diff.pTrans = TRUE,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          diff.pExit=TRUE,
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 12)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),3)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")$state,c("A","C","A"))

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.when.infected = getHostData(test.nosoiA, "table.hosts") %>% group_by(inf.in) %>% summarise(N=length(hosts.ID))
  expect_equal(subset(Where.when.infected, inf.in == "B")$N,integer(0))
  expect_equal(subset(Where.when.infected, inf.in == "A")$N,1)
  expect_equal(subset(Where.when.infected, inf.in == "C")$N,968)

  Where.when.died = subset(getHostData(test.nosoiA, "table.hosts"),active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  expect_equal(subset(Where.when.died, current.in == "C")$N,integer(0))
  expect_equal(subset(Where.when.died, current.in == "D")$N,207)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})

test_that("Epidemic dying out", {
  skip_if_not_installed("igraph")
  library(igraph)

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

  set.seed(10)
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

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),1)
  expect_equal(test.nosoiA$total.time,4)

  skip_if_not_installed("dplyr")
  skip_if_not_installed("magrittr")
  library(dplyr)
  library(magrittr)

  Where.at.end = getHostData(test.nosoiA, "table.hosts") %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.at.end, current.in == "A")$N,1)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})


test_that("Movement is coherent with single introduction, no pMove, no die, diff nContact with hostCount needed", {
  skip_if_not_installed("igraph")
  library(igraph)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  p_Exit_fct  <- function(t){return(0.05)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact <- function(t, current.in, host.count){

    temp.val = 30 - host.count

    if(temp.val <= 0) {
      return(0)
    }
    if(temp.val >= 0) {
      if(current.in=="A"){
        return(round((temp.val/30)*rnorm(1, 3, 1), 0))}
      if(current.in=="B"){return(0)}
      if(current.in=="C"){
        return(round((temp.val/30)*rnorm(1, 6, 1), 0))}
    }
  }

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(1050)
  test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                          length=100,
                          max.infected=200,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.nContact=TRUE,
                          hostCount.nContact=TRUE,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(getHostData(test.nosoiA, "table.hosts")[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 10)

  expect_equal(nrow(getHostData(test.nosoiA, "table.hosts")),204)

  #Movement
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-2")),5)
  expect_equal(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-2")$state,c("A","C","A","C","B"))

  #Number of host at each loc at each time
  out = getDynamic(test.nosoiA) #  ggplot(out, aes(x=t,y=Count,color=state)) + geom_line() + geom_hline(yintercept=30)
  expect_equal(subset(out,state=="C" & t==26)$Count,28)
  expect_equal(subset(out,state=="B" & t==26)$Count,15)
  expect_equal(subset(out,state=="A" & t==25)$Count,30)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0$N.inactive,
               ifelse(length(r_0$R0.dist) == 1 && is.na(r_0$R0.dist), 0, length(r_0$R0.dist)))
})
