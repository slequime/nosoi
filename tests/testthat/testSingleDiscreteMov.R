context("Checking the movement in discrete structure")

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
    test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                            length=20,
                            max.infected=100,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            timeContact=time_contact,
                            param.timeContact=NA,
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
    test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                            length=20,
                            max.infected=100,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            timeContact=time_contact,
                            param.timeContact=NA,
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
    test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                            length=20,
                            max.infected=100,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            timeContact=time_contact,
                            param.timeContact=NA,
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
    test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                            length=20,
                            max.infected=100,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            timeContact=time_contact,
                            param.timeContact=NA,
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
    test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                            length=20,
                            max.infected=100,
                            init.individuals=1,
                            init.structure="D",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            timeContact=time_contact,
                            param.timeContact=NA,
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
test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                        length=20,
                        max.infected=100,
                        init.individuals=1,
                        init.structure="A",
                        structure.matrix=transition.matrix,
                        pMove=p_Move_fct,
                        param.pMove=NA,
                        timeContact=time_contact,
                        param.timeContact=NA,
                        pTrans = proba,
                        param.pTrans = list(p_max=p_max_fct,
                                            t_incub=t_incub_fct),
                        pExit=p_Exit_fct,
                        param.pExit=NA
)

#Structure
g <- graph.data.frame(test.nosoiA$table.hosts[,c(1,2)],directed=F)
expect_equal(transitivity(g, type="global"), 0)
expect_equal(clusters(g, "weak")$no, 1)
expect_equal(diameter(g, directed=F, weights=NA), 6)

#Movement
expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-3")),3)
expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-3")$state,c("A","C","B"))

Where.at.end = test.nosoiA$table.hosts %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

expect_equal(subset(Where.at.end, current.in == "A")$N,62)
expect_equal(subset(Where.at.end, current.in == "B")$N,17)
expect_equal(subset(Where.at.end, current.in == "C")$N,43)
})

#--------------------------------------------------------------------------------------------------------------------------------------------

test_that("Movement is coherent with single introduction, complex pMove", {
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
test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                        length=20,
                        max.infected=1000,
                        init.individuals=1,
                        init.structure="A",
                        structure.matrix=transition.matrix,
                        pMove=p_Move_fct,
                        param.pMove=list(pMove.param1=p_Move_param1_fct),
                        timeContact=time_contact,
                        param.timeContact=NA,
                        pTrans = proba,
                        param.pTrans = list(p_max=p_max_fct,
                                            t_incub=t_incub_fct),
                        pExit=p_Exit_fct,
                        param.pExit=NA
)

#Structure
g <- graph.data.frame(test.nosoiA$table.hosts[,c(1,2)],directed=F)
expect_equal(transitivity(g, type="global"), 0)
expect_equal(clusters(g, "weak")$no, 1)
expect_equal(diameter(g, directed=F, weights=NA), 6)

#Movement
expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-1")),11)
expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-1")$state,c("A","C","A","C","A","B","C","B","C","A","C"))

Where.at.end = test.nosoiA$table.hosts %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

expect_equal(subset(Where.at.end, current.in == "A")$N,83)
expect_equal(subset(Where.at.end, current.in == "B")$N,19)
expect_equal(subset(Where.at.end, current.in == "C")$N,16)
})

test_that("Movement is coherent with single introduction, constant but different pMove, 1 loc (C) is sink. Ce tombeau sera votre tombeau !", {

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
test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                        length=20,
                        max.infected=10000,
                        init.individuals=1,
                        init.structure="A",
                        structure.matrix=transition.matrix,
                        diff.pMove=TRUE,
                        pMove=p_Move_fct,
                        param.pMove=param.pMove,
                        timeContact=time_contact,
                        param.timeContact=NA,
                        pTrans = proba,
                        param.pTrans = list(p_max=p_max_fct,
                                            t_incub=t_incub_fct),
                        pExit=p_Exit_fct,
                        param.pExit=NA
)

#Structure
g <- graph.data.frame(test.nosoiA$table.hosts[,c(1,2)],directed=F)
expect_equal(transitivity(g, type="global"), 0)
expect_equal(clusters(g, "weak")$no, 1)
expect_equal(diameter(g, directed=F, weights=NA), 6)

#Movement
expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-1")),2)
expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-1")$state,c("A","C"))

Where.at.end = test.nosoiA$table.hosts %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
expect_equal(subset(Where.at.end, current.in == "C")$N,138)
})

test_that("Error message pops up if different pMove poorly formated", {
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
  test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                          length=20,
                          max.infected=10000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          diff.pMove=TRUE,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          timeContact=time_contact,
                          param.timeContact=NA,
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
  test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                          length=20,
                          max.infected=10000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          diff.pMove=TRUE,
                          pMove=p_Move_fct,
                          param.pMove=param.pMove,
                          timeContact=time_contact,
                          param.timeContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(test.nosoiA$table.hosts[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 6)

  #Movement
  expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-1")),2)
  expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-1")$state,c("A","C"))

  Where.at.end = test.nosoiA$table.hosts %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  expect_equal(subset(Where.at.end, current.in == "C")$N,52)
})
