context("Parsing functions")

test_that("Missing args and error messages", {

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t, p_max, t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  p_exit_param1 <- function(x){rnorm(x,mean = 10,sd=2)}

  p_Exit_fct_complex  <- function(t, pExit.param1){plogis(t, pExit.param1, 2)}

  expect_error(
    test.nosoiA <- nosoiSim(type = "single", popStructure="none",
                            length = 40,
                            max.infected = 100,
                            init.individuals = 1,
                            nContact = time_contact,
                            pTrans = proba,
                            pExit = p_Exit_fct,
                            param.pExit = NA
    ),
    "argument \"param.pTrans\" is missing, with no default"
  )

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="none",
                            length=40,
                            max.infected=100,
                            init.individuals=1,
                            nContact=time_contact,
                            pTrans = proba,
                            param.pTrans = NA,
                            pExit = p_Exit_fct,
                            param.pExit = NA
    ),
    "There is a probleme with your function pTrans: you should provide a parameter list named param.pTrans."
  )

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="none",
                            length=40,
                            max.infected=100,
                            init.individuals=1,
                            nContact=time_contact,
                            pTrans = proba,
                            param.pTrans = list(p_max = p_max_fct,
                                                t_incub = t_incub_fct),
                            pExit = p_Exit_fct
    ),
    "argument \"param.pExit\" is missing, with no default"
  )

  expect_error(
    test.nosoiA <- nosoiSim(type="single", popStructure="none",
                            length=40,
                            max.infected=100,
                            init.individuals=1,
                            nContact=time_contact,
                            pTrans = proba,
                            param.pTrans = list(p_max = p_max_fct,
                                                t_incub = t_incub_fct),
                            pExit = p_Exit_fct_complex,
                            param.pExit = NA
    ),
    "There is a probleme with your function pExit: you should provide a parameter list named param.pExit."
  )

  proba2 <- 0.5
  expect_error(
    nosoiSim(type="single", popStructure="none",
                            length=40,
                            max.infected=100,
                            init.individuals=1,
                            nContact=time_contact,
                            pTrans = proba2,
                            param.pTrans = list(p_max = p_max_fct,
                                                t_incub = t_incub_fct),
                            pExit = p_Exit_fct,
                            param.pExit = NA
    ),
    "You must specify pTrans as a function.", fixed=TRUE
  )

  proba3 <- function(z){return(0.08)}
  expect_error(
    nosoiSim(type="single", popStructure="none",
             length=40,
             max.infected=100,
             init.individuals=1,
             nContact=time_contact,
             pTrans = proba3,
             param.pTrans = list(p_max = p_max_fct,
                                 t_incub = t_incub_fct),
             pExit = p_Exit_fct,
             param.pExit = NA
    ),
    "pTrans must be a function of 't', placed as a first argument of the function.", fixed=TRUE
  )

  expect_error(
    nosoiSim(type="single", popStructure="none",
             length=40,
             max.infected=100,
             init.individuals=1,
             nContact=time_contact,
             pTrans = proba,
             param.pTrans = list(p_max2 = p_max_fct,
                                 t_incub = t_incub_fct),
             pExit = p_Exit_fct,
             param.pExit = NA
    ),
    "Parameter name in param.pTrans should match the name used in pTrans.", fixed=TRUE
  )

  expect_error(
    nosoiSim(type="single", popStructure="none",
             length=40,
             max.infected=100,
             init.individuals=1,
             nContact=time_contact,
             pTrans = proba,
             timeDep.pTrans=TRUE,
             param.pTrans = list(p_max = p_max_fct,
                                 t_incub = t_incub_fct),
             pExit = p_Exit_fct,
             param.pExit = NA
    ),
    "pTrans should have 'prestime' as the second variable. timeDep.pTrans is TRUE.", fixed=TRUE
  )

})

test_that("Error message with discrete structure", {

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}

   p_Move_param1_fct <- function(x){rnorm(x,mean = 10,sd=2)}

   p_Move_fct_ <- function(t,pMove.param1,current.in){
     if(current.in=="A"){return(plogis(t,1+pMove.param1,2))}
     if(current.in=="B"){return(plogis(t,2+pMove.param1,2))}
     if(current.in=="C"){return(plogis(t,pMove.param1,2)/1000)}}

   p_Move_fct_uncorrect  <- function(t,pMove.param1,current.in){
     if(current.in=="A"){return(plogis(t,1+pMove.param1,2))}
     if(current.in=="B"){return(plogis(t,2+pMove.param1,2))}
     if(current.in=="C"){return(plogis(t,pMove.param1,2)/1000)}}

   p_Move_fct_uncorrect2  <- function(t,prestime,pMove.param1,current.in){
     if(current.in=="A"){return(plogis(prestime,1+pMove.param1,2))}
     if(current.in=="B"){return(plogis(prestime,2+pMove.param1,2))}
     if(current.in=="C"){return(plogis(prestime,pMove.param1,2)/1000)}}

   param.pMove = list(pMove.param1=p_Move_param1_fct)

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  transition.matrix = matrix(c(0,0.2,0.4,0.5,0,0.6,0.5,0.8,0),nrow = 3, ncol = 3,dimnames=list(c("A","B","C"),c("A","B","C")))

  set.seed(750)
  expect_error(
    nosoiSim(type="single", popStructure="discrete",
             length=20,
             max.infected=10000,
             init.individuals=1,
             init.structure="A",
             structure.matrix=transition.matrix,
             diff.pMove=TRUE,
             timeDep.pMove=TRUE,
             pMove=p_Move_fct_uncorrect2,
             param.pMove=param.pMove,
             nContact=time_contact,
             param.nContact=NA,
             pTrans = proba,
             param.pTrans = list(p_max=p_max_fct,
                                 t_incub=t_incub_fct),
             pExit=p_Exit_fct,
             param.pExit=NA
    ),"pMove should have 'current.in' as the third variable. diff.pMove is TRUE."
  )

  expect_error(
    nosoiSim(type="single", popStructure="discrete",
             length=20,
             max.infected=10000,
             init.individuals=1,
             init.structure="A",
             structure.matrix=transition.matrix,
             diff.pMove=TRUE,
             pMove=p_Move_fct_uncorrect2,
             param.pMove=param.pMove,
             nContact=time_contact,
             param.nContact=NA,
             pTrans = proba,
             param.pTrans = list(p_max=p_max_fct,
                                 t_incub=t_incub_fct),
             pExit=p_Exit_fct,
             param.pExit=NA
    ),"pMove should have 'current.in' as the second variable. diff.pMove is TRUE."
  )
})

test_that("Error message with continuous structure", {

  library(raster)

  #Generating a raster for the movement
  set.seed(10)

  test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
  test.raster[] <- runif(10000, -80, 180)
  test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
  # plot(test.raster)

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_Move_fct  <- function(t){return(0.1)}

  sd_Move_param1_fct <- function(x){rnorm(x,mean = 10,sd=2)}

  sdMove_fct = function(t,prestime,sdMove.param1,current.env.value){return(100/(current.env.value+1))}

  param.sdMove = list(sdMove.param1=sd_Move_param1_fct)

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=0.3}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  start.pos <- c(0,0)

  expect_error(
    nosoiSim(type="single", popStructure="continuous",
                          length=200,
                          max.infected=500,
                          init.individuals=1,
                          init.structure=start.pos,
                          structure.raster=test.raster,
                          pMove=p_Move_fct,
                          param.pMove=NA,
                          diff.sdMove=TRUE,
                          timeDep.sdMove=TRUE,
                          sdMove=sdMove_fct,
                          param.sdMove=param.sdMove,
                          attracted.by.raster=FALSE,
                          nContact=time_contact,
                          param.nContact=NA,
                          pTrans = proba,
                          param.pTrans = list(t_incub=t_incub_fct),
                          pExit=p_Exit_fct,
                          param.pExit=NA
  ),"sdMove should have 'current.env.value' as the third variable. diff.sdMove is TRUE."
)


})


test_that("Actual parsing", {

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t, p_max, t_incub){
    if(t <= t_incub){p = 0}
    if(t >= t_incub){p = p_max}
    return(p)
  }

  parsedProba <- parseFunction(proba, list(p_max = p_max_fct, t_incub = t_incub_fct), "blabla")
  expect_equal(parsedProba$vectArgs, c("inf.time", "p_max", "t_incub"))

  parsedProba <- parseFunction(p_Exit_fct, NA, "blabla")
  expect_equal(parsedProba$vectArgs, "inf.time")

})


