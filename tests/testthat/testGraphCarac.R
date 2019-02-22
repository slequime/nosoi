context("Graph of the transmission tree structure")

test_that("Transmission is coherent with single introduction", {
  library(igraph)
  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(x){0.08}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(x){round(rnorm(x, 3, 1), 0)}
  #Test with one introduction
  set.seed(805)
  test.nosoiA <- nosoiSim(type="single",geo="none",
                         length=40,
                         max.infected=100,
                         init.individuals=1,
                         timeContact=time_contact,
                         pTrans = proba,
                         param.pTrans = list(p_max=p_max_fct,
                                             t_incub=t_incub_fct),
                         pExit=p_Exit_fct,
                         param.pExit=NA
                         )
  g <- graph.data.frame(test.nosoiA[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 1)
  expect_equal(diameter(g, directed=F, weights=NA), 5)
})

test_that("Transmission is coherent with multiple introduction", {
  library(igraph)
  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(x){0.08}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(x){round(rnorm(x, 3, 1), 0)}
  #Test with one introduction
  set.seed(110)
  test.nosoiB <- nosoiSim(type="single",geo="none",
                          length=40,
                          max.infected=100,
                          init.individuals=3,
                          timeContact=time_contact,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  )
  g <- graph.data.frame(test.nosoiB[,c(1,2)],directed=F)

  expect_equal(transitivity(g, type="global"), 0)
  expect_equal(clusters(g, "weak")$no, 3)
  expect_equal(diameter(g, directed=F, weights=NA), 4)
})
