context("Checking the differential probabilities in discrete structure")

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
  test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                          length=20,
                          max.infected=1000,
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
                          diff.pExit=TRUE,
                          pExit=p_Exit_fct,
                          param.pExit=NA
  ),
  "pExit should have a realisation for each possible state. diff.pExit == TRUE."
  )

  p_Exit_fct  <- function(x){return(0.08)}

  time_contact <- function(t,current.in){
    if(current.in=="A"){
      return(round(rnorm(1, 3, 1), 0))}
    if(current.in=="C"){
      return(round(rnorm(1, 6, 1), 0))}
  }

  expect_error(
    test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                            length=20,
                            max.infected=1000,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            diff.timeContact=TRUE,
                            timeContact=time_contact,
                            param.timeContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=NA,
                            pExit=p_Exit_fct,
                            param.pExit=NA
    ),
    "timeContact should have a realisation for each possible state. diff.timeContact == TRUE."
  )

  time_contact = function(x){round(rnorm(x, 3, 1), 0)}

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
    test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                            length=20,
                            max.infected=1000,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            timeContact=time_contact,
                            param.timeContact=NA,
                            diff.pTrans=TRUE,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            pExit=p_Exit_fct,
                            param.pExit=NA
    ),
    "pTrans should have a realisation for each possible state. diff.pTrans == TRUE."
  )

})


test_that("Movement is coherent with single introduction, constant pMove, diff pExit", {
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
test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                        length=20,
                        max.infected=1000,
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
                        diff.pExit=TRUE,
                        pExit=p_Exit_fct,
                        param.pExit=NA
)

#Structure
g <- graph.data.frame(test.nosoiA$table.hosts[,c(1,2)],directed=F)
expect_equal(transitivity(g, type="global"), 0)
expect_equal(clusters(g, "weak")$no, 1)
expect_equal(diameter(g, directed=F, weights=NA), 7)

#Movement
expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-3")),2)
expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-3")$state,c("A","C"))

Where.when.exit = subset(test.nosoiA$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))

expect_equal(subset(Where.when.exit, current.in == "A")$N,integer(0))
expect_equal(subset(Where.when.exit, current.in == "B")$N,44)
expect_equal(subset(Where.when.exit, current.in == "C")$N,34)
})

test_that("Movement is coherent with single introduction, constant pMove, diff pTrans ", {
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
  test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                          length=20,
                          max.infected=1000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          timeContact=time_contact,
                          param.timeContact=NA,
                          diff.pTrans = TRUE,
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
  expect_equal(diameter(g, directed=F, weights=NA), 7)

  #Movement
  expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-3")),2)
  expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-3")$state,c("A","C"))

  Where.when.infected = subset(test.nosoiA$table.hosts,active==0) %>% group_by(inf.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.infected, inf.in == "B")$N,integer(0))
  expect_equal(subset(Where.when.infected, inf.in == "A")$N,56)
  expect_equal(subset(Where.when.infected, inf.in == "C")$N,17)
})

test_that("Movement is coherent with single introduction, constant pMove, diff timeContact ", {
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
  test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                          length=20,
                          max.infected=1000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.timeContact=TRUE,
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
  expect_equal(diameter(g, directed=F, weights=NA), 8)

  #Movement
  expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-3")),2)
  expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-3")$state,c("A","C"))

  Where.when.infected = subset(test.nosoiA$table.hosts,active==0) %>% group_by(inf.in) %>% summarise(N=length(hosts.ID))

  expect_equal(subset(Where.when.infected, inf.in == "B")$N,integer(0))
  expect_equal(subset(Where.when.infected, inf.in == "A")$N,27)
  expect_equal(subset(Where.when.infected, inf.in == "C")$N,67)
})

test_that("Movement is coherent with single introduction, all parameters are diff", {
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
  test.nosoiA <- nosoiSim(type="single",structure=TRUE,
                          length=40,
                          max.infected=1000,
                          init.individuals=1,
                          init.structure="A",
                          structure.matrix=transition.matrix,
                          diff.pMove=TRUE,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.timeContact=TRUE,
                          timeContact=time_contact,
                          param.timeContact=NA,
                          diff.pTrans = TRUE,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          diff.pExit=TRUE,
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )

  #Structure
  g <- graph.data.frame(test.nosoiA$table.hosts[,c(1,2)],directed=F)
  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 12)

  #Movement
  expect_equal(nrow(subset(test.nosoiA$table.state, hosts.ID == "H-1")),3)
  expect_equal(subset(test.nosoiA$table.state, hosts.ID == "H-1")$state,c("A","C","A"))

  Where.when.infected = test.nosoiA$table.hosts %>% group_by(inf.in) %>% summarise(N=length(hosts.ID))
  expect_equal(subset(Where.when.infected, inf.in == "B")$N,integer(0))
  expect_equal(subset(Where.when.infected, inf.in == "A")$N,1)
  expect_equal(subset(Where.when.infected, inf.in == "C")$N,968)

  Where.when.died = subset(test.nosoiA$table.hosts,active==0) %>% group_by(current.in) %>% summarise(N=length(hosts.ID))
  expect_equal(subset(Where.when.died, current.in == "C")$N,integer(0))
  expect_equal(subset(Where.when.died, current.in == "D")$N,207)
})
