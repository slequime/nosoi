context("Checking movement and differential expressions in continuous space")

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

  moveDist_fct = function(t){return(1)}

  p_Exit_fct  <- function(t){return(0.08)}

  proba <- function(t,p_max,t_incub){
    if(t <= t_incub){p=0}
    if(t >= t_incub){p=p_max}
    return(p)
  }

  time_contact = function(t){round(rnorm(1, 3, 1), 0)}

  start.pos <- c(0,0)

  expect_error(
    test.nosoiA <- nosoiSim(type="single", structure=TRUE, continuous = TRUE,
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            moveDist=moveDist_fct,
                            param.moveDist=NA,
                            attracted.by.raster=TRUE,
                            timeContact=time_contact,
                            param.timeContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "pExit should have 'current.env.value' as a variable. diff.pExit == TRUE."
  )

  p_Exit_fct  <- function(t,current.env.value){(1-current.env.value)}

  expect_error(
    test.nosoiA <- nosoiSim(type="single", structure=TRUE, continuous = TRUE,
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            diff.pMove=TRUE,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            moveDist=moveDist_fct,
                            param.moveDist=NA,
                            attracted.by.raster=TRUE,
                            timeContact=time_contact,
                            param.timeContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "pMove should have 'current.env.value' as a variable. diff.pMove == TRUE."
  )

  p_Move_fct  <- function(t,current.env.value){current.env.value/1000}

  expect_error(
    test.nosoiA <- nosoiSim(type="single", structure=TRUE, continuous = TRUE,
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            diff.pMove=TRUE,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            diff.moveDist=TRUE,
                            moveDist=moveDist_fct,
                            param.moveDist=NA,
                            attracted.by.raster=TRUE,
                            timeContact=time_contact,
                            param.timeContact=NA,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "moveDist should have 'current.env.value' as a variable. diff.moveDist == TRUE."
  )

  moveDist_fct = function(t,current.env.value){return(100/current.env.value+1)}

  expect_error(
    test.nosoiA <- nosoiSim(type="single", structure=TRUE, continuous = TRUE,
                            length=365,
                            max.infected=10000,
                            init.individuals=1,
                            init.structure=start.pos,
                            structure.raster=test.raster,
                            diff.pMove=TRUE,
                            pMove=p_Move_fct,
                            param.pMove=NA,
                            diff.moveDist=TRUE,
                            moveDist=moveDist_fct,
                            param.moveDist=NA,
                            attracted.by.raster=TRUE,
                            timeContact=time_contact,
                            param.timeContact=NA,
                            diff.pTrans=TRUE,
                            pTrans = proba,
                            param.pTrans = list(p_max=p_max_fct,
                                                t_incub=t_incub_fct),
                            diff.pExit=TRUE,
                            pExit=p_Exit_fct,
                            param.pExit=NA),
    "pTrans should have 'current.env.value' as a variable. diff.pTrans == TRUE."
  )
})





































#
#
#
#
#
#
#
#
#
#
# test_that("Movement is coherent with single introduction, constant pMove", {
#
#   #Generating a raster the for movement
#   set.seed(860)
#
#   test.raster <- raster(nrows=100, ncols=100, xmn=-50, xmx=50, ymn=-50,ymx=50)
#   test.raster[] <- runif(10000, -80, 180)
#   test.raster <- focal(focal(test.raster, w=matrix(1, 5, 5), mean), w=matrix(1, 5, 5), mean)
#   # plot(test.raster)
#
#   t_incub_fct <- function(x){rnorm(x,mean = 5,sd=1)}
#   p_max_fct <- function(x){rbeta(x,shape1 = 5,shape2=2)}
#   p_Move_fct  <- function(t){return(0.1)}
#
#   moveDist_fct = function(t,current.env.value){return(100/(current.env.value+1))}
#
#   p_Exit_fct  <- function(t){return(0.08)}
#
#   proba <- function(t,p_max,t_incub){
#     if(t <= t_incub){p=0}
#     if(t >= t_incub){p=p_max}
#     return(p)
#   }
#
#   time_contact = function(t){round(rnorm(1, 3, 1), 0)}
#
#   start.pos <- c(0,0)
#
#   test.nosoiA <- nosoiSim(type="single", structure=TRUE, continuous = TRUE,
#                           length=365,
#                           max.infected=10000,
#                           init.individuals=1,
#                           init.structure=start.pos,
#                           structure.raster=test.raster,
#                           pMove=p_Move_fct,
#                           param.pMove=NA,
#                           diff.moveDist=TRUE,
#                           moveDist=moveDist_fct,
#                           param.moveDist=NA,
#                           attracted.by.raster=TRUE,
#                           timeContact=time_contact,
#                           param.timeContact=NA,
#                           pTrans = proba,
#                           param.pTrans = list(p_max=p_max_fct,
#                                               t_incub=t_incub_fct),
#                           pExit=p_Exit_fct,
#                           param.pExit=NA)
#
# ggplot() +
#     layer_spatial(tc) +
#     scale_fill_gradient(low = "gray95", high = "forestgreen", limits=c(0,NA),na.value="white", name=NULL) +
#     # scale_fill_gradient(low = "cornsilk", high = "darkgoldenrod4",na.value=NA, name=NULL) +
#     # scale_y_continuous(breaks = seq(30, 90, by = 5), labels = NULL) +
#     theme(panel.background = element_blank(),
#           panel.ontop = TRUE,
#           panel.grid.major = element_line(size = 0.25, linetype = 'dashed',colour = "grey70"),
#           axis.ticks=element_blank(),
#           plot.caption=element_text(size=10)) +
#     labs(caption = "Elevation (source: Spatial Data Access Tool - SDAT)") +
#     geom_spatial_point(data=test.nosoiA$table.state,aes(state.x,state.y),color="firebrick3",alpha=0.1)
#
#
# })
