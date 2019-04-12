context("Parsing functions")

test_that("Missing args", {

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(x){rep(0.08,length(x))}

  proba <- function(t, p_max, t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(x){round(rnorm(x, 3, 1), 0)}

  p_exit_param1 <- function(x){rnorm(x,mean = 10,sd=2)}

  p_Exit_fct_complex  <- function(t, pExit.param1){plogis(t, pExit.param1, 2)}

  expect_error(
    test.nosoiA <- nosoiSim(type = "single", structure = FALSE,
                            length = 40,
                            max.infected = 100,
                            init.individuals = 1,
                            timeContact = time_contact,
                            pTrans = proba,
                            pExit = p_Exit_fct,
                            param.pExit = NA
    ),
    "argument \"param.pTrans\" is missing, with no default"
  )

  expect_error(
    test.nosoiA <- nosoiSim(type="single",structure=FALSE,
                            length=40,
                            max.infected=100,
                            init.individuals=1,
                            timeContact=time_contact,
                            pTrans = proba,
                            param.pTrans = NA,
                            pExit = p_Exit_fct,
                            param.pExit = NA
    ),
    "There is a probleme with your function pTrans: you should provide a parameter list named param.pTrans."
  )

  expect_error(
    test.nosoiA <- nosoiSim(type="single",structure=FALSE,
                            length=40,
                            max.infected=100,
                            init.individuals=1,
                            timeContact=time_contact,
                            pTrans = proba,
                            param.pTrans = list(p_max = p_max_fct,
                                                t_incub = t_incub_fct),
                            pExit = p_Exit_fct
    ),
    "argument \"param.pExit\" is missing, with no default"
  )

  expect_error(
    test.nosoiA <- nosoiSim(type="single",structure=FALSE,
                            length=40,
                            max.infected=100,
                            init.individuals=1,
                            timeContact=time_contact,
                            pTrans = proba,
                            param.pTrans = list(p_max = p_max_fct,
                                                t_incub = t_incub_fct),
                            pExit = p_Exit_fct_complex,
                            param.pExit = NA
    ),
    "There is a probleme with your function pExit: you should provide a parameter list named param.pExit."
  )

})

test_that("Actual parsing", {

  t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
  p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
  p_Exit_fct  <- function(x){rep(0.08,length(x))}

  proba <- function(t, p_max, t_incub){
    if(t <= t_incub){p = 0}
    if(t >= t_incub){p = p_max}
    return(p)
  }

  parsedProba <- parseFunction(proba, list(p_max = p_max_fct, t_incub = t_incub_fct), "blabla")
  expect_equal(parsedProba$type, "complex")
  expect_equal(parsedProba$nArgs, 2)
  expect_equal(parsedProba$vectArgs, c("inf.time", "p_max", "t_incub"))

  parsedProba <- parseFunction(p_Exit_fct, NA, "blabla")
  expect_equal(parsedProba$type, "simple")
  expect_equal(parsedProba$nArgs, 0)
  expect_equal(parsedProba$vectArgs, "inf.time")

})

# test_that("Actual parsing with diff", {
#
#   p_Move_fct  <- function(x,state){
#     if(state=="A"){return(0.1)}
#     if(state=="B"){return(0.3)}
#     if(state=="C"){return(0.01)}}
#
#   param.pMove = NA
#
#   parsedProba <- parseFunction(p_Move_fct, param.pMove, "blabla", diff=TRUE)
#   expect_equal(parsedProba$type, "simple-Diff")
#   expect_equal(parsedProba$nArgs, 0)
#   expect_equal(parsedProba$vectArgs, c("inf.time","state"))
# })
