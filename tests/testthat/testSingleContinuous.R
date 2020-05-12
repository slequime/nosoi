context("Testing single-host with continuous structure")

test_that("Error message pops out when missing state in diff functions", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)

  # plot(test.raster)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  sdMove_fct = function(t){return(1)}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  start.pos <- c(0,0)

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            sdMove=sdMove_fct,
                            param.sdMove=NA,
                            attracted.by.raster=TRUE,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "Your are missing some function argument in pExit. diff and/or timeDep.pExit is/are TRUE."
  )

  p_Exit_fct  <- function(t,current.env.value){(1-current.env.value)}

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            diff.pMove=TRUE,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            sdMove=sdMove_fct,
                            param.sdMove=NA,
                            attracted.by.raster=TRUE,
                            nContact=time_contact,
                            param.nContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "Your are missing some function argument in pMove. diff and/or timeDep.pMove is/are TRUE."
  )

  p_Move_fct  <- function(t,current.env.value){current.env.value/1000}

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            diff.pMove=TRUE,
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
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "Your are missing some function argument in sdMove. diff and/or timeDep.sdMove is/are TRUE."
  )

  sdMove_fct = function(t,current.env.value){return(100/current.env.value+1)}

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            diff.pMove=TRUE,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            diff.sdMove=TRUE,
                            sdMove=sdMove_fct,
                            param.sdMove=NA,
                            attracted.by.raster=TRUE,
                            nContact=time_contact,
                            param.nContact=NA,
                            diff.pTrans=TRUE,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "pTrans should have 'current.env.value' as the second variable. diff.pTrans is TRUE."
  )
})

test_that("Diffusion in continuous space", {
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
  # plot(test.raster)

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

  test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                          length=200,
                          max.infected=500,
                          init.individuals=1,
                          init.structure=start.pos,
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
                          param.pExit=NA)

  expect_equal(nrow(getHostData(test.nosoiA, "table.hosts")),648)
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),3)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0_old <- getR0Old(test.nosoiA)
  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0_old$R0.mean, r_0$R0.mean)
})

test_that("Epidemic dying out", {
  library(raster)

  #Generating a raster the for movement
  set.seed(10)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
  # plot(test.raster)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_Move_fct  <- function(t){return(0.1)}

  sdMove_fct = function(t,current.env.value){return(100/(current.env.value+1))}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=0.3}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  start.pos <- c(0,0)

  test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                          length=200,
                          max.infected=500,
                          init.individuals=1,
                          init.structure=start.pos,
                          structure.raster=test.raster,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.sdMove=TRUE,
                          sdMove=sdMove_fct,
                          param.sdMove=NA,
                          attracted.by.raster=FALSE,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA)

  expect_equal(nrow(getHostData(test.nosoiA, "table.hosts")),3)
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),3)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0_old <- getR0Old(test.nosoiA)
  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0_old$R0.mean, r_0$R0.mean)
})

test_that("Diffusion in continuous space with host count", {
  # library(profvis)

  # profvis({
  library(raster)

  #Generating a raster the for movement
  set.seed(860)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
  # plot(test.raster)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Move_fct  <- function(t){return(0.1)}

  sdMove_fct = function(t){return(0.2)}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact <- function(t, current.env.value, host.count){

    temp.val = round(((current.env.value-host.count)/current.env.value)*rnorm(1, 3, 1), 0)

    if(length(temp.val) == 0 || temp.val <= 0) {
      return(0)
    }
    if(temp.val >= 0) {
      return(temp.val)
    }
  }

  start.pos <- c(0,0)

  set.seed(9897)
  test.nosoiA <- nosoiSim(type="single", popStructure="continuous",
                          length=20,
                          max.infected=100,
                          init.individuals=1,
                          init.structure=start.pos,
                          structure.raster=test.raster,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.sdMove=FALSE,
                          sdMove=sdMove_fct,
                          param.sdMove=NA,
                          attracted.by.raster=FALSE,
                          nContact=time_contact,
                          diff.nContact=TRUE,
                          hostCount.nContact=TRUE,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(p_max=p_max_fct,
                                              t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA)

  expect_equal(nrow(getHostData(test.nosoiA, "table.hosts")),115)
  expect_equal(nrow(subset(getHostData(test.nosoiA, "table.state"), hosts.ID == "H-1")),3)

  skip_if_not_installed("dplyr")
  dynOld <- getDynamicOld(test.nosoiA)
  dynNew <- getDynamic(test.nosoiA)
  expect_equal(dynOld, dynNew)

  r_0_old <- getR0Old(test.nosoiA)
  r_0 <- getR0(test.nosoiA)
  expect_equal(r_0_old$R0.mean, r_0$R0.mean)
})
