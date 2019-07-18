context("Other error messages")

test_that("Error messages on NosoiSim", {

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

  expect_error(
  test.nosoiC <- nosoiSim(type="triple", popStructure="none",
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
  ),
  "Type of transmission should be 'single' or 'dual'-host."
  )

  expect_error(
    test.nosoiC <- nosoiSim(type="single", popStructure="multiple",
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
    ),
    "Unrecognized parameters for population structure, should be 'none,'discrete' or 'continuous'."
  )

})

test_that("Error messages on SanityCheck", {
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

  expect_error(
    nosoiSim(type="single", popStructure="none",
                            length=1,
                            max.infected=100,
                            init.individuals=3,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            pExit=p_Exit_fct,
                            param.pExit = list(pExit.param1=p_exit_param1)
    ),
    "You must specify a length (in time units) for your simulation (bigger than 1).",fixed=TRUE
  )

  expect_error(
    test.nosoiC <- nosoiSim(type="single", popStructure="none",
                            length=100,
                            max.infected=NA,
                            init.individuals=3,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            pExit=p_Exit_fct,
                            param.pExit = list(pExit.param1=p_exit_param1)
    ),
    "You must specify a maximum number of infected hosts (bigger than 1).",fixed=TRUE
  )

  expect_error(
    test.nosoiC <- nosoiSim(type="single",  popStructure="none",
                            length=100,
                            max.infected=100,
                            init.individuals=NA,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            pExit=p_Exit_fct,
                            param.pExit = list(pExit.param1=p_exit_param1)
    ),
    "The transmission chain should be started by 1 or more (integer) individuals.",fixed=TRUE
  )

})

test_that("Error messages on MatrixSanityCheck", {

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
  expect_error(
    test.nosoiA <- nosoiSim(type="dual", popStructure="discrete",
                          length.sim=20,
                          max.infected.A=1000,
                          max.infected.B=1000,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A="D",
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
                          prefix.host.B="V"),
  "init.structure should be a state present in structure.matrix.")


})

test_that("Error messages on RasterSanityCheck", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
  # plot(test.raster)
  bad.raster <- "raster"

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
  bad.start.pos <- c(100,100)

  expect_error(nosoiSim(type="single", popStructure="continuous",
                          length=200,
                          max.infected=500,
                          init.individuals=1,
                          init.structure=start.pos,
                          structure.raster=bad.raster,
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
                          param.pExit=NA),
  "structure.raster must be a raster (class RasterLayer).", fixed=TRUE)

  expect_error(nosoiSim(type="single", popStructure="continuous",
                        length=200,
                        max.infected=500,
                        init.individuals=1,
                        init.structure=bad.start.pos,
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
                        param.pExit=NA),
               "Your starting position (init.structure) should be on the raster.", fixed=TRUE)


  set.seed(805)
  expect_error(
    nosoiSim(type="dual", popStructure="continuous",
                          length.sim=200,
                          max.infected.A=500,
                          max.infected.B=500,
                          init.individuals.A=1,
                          init.individuals.B=0,
                          init.structure.A=bad.start.pos,
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
                          prefix.host.B="V"),
    "Your starting position (init.structure) should be on the raster.",fixed=TRUE)

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="discrete",
                            length=100,
                            max.infected=200,
                            init.individuals=1,
                            init.structure="A",
                            structure.matrix=transition.matrix,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            diff.nContact=FALSE,
                            hostCount.nContact=TRUE,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            pExit=p_Exit_fct,
                            param.pExit=NA
    ),
    "diff.nContact should be TRUE to use hostCount.nContact.")
  })

test_that("Error messages on drawBernouilli", {
  test = data.frame(A=0.5,B=0.6,C=0.1)
  expect_error(drawBernouilli(test),
               "Function 'drawBernouilli' should be applied to a vector.")
})
