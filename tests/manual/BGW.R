context("Simple case equivalent to the BGW process")

test_that("Distribution of the number of active individuals - surcritic", {
  ## Parameters
  n_pop <- 100000                # Population size
  p_contact <- 0.001             # Contact proba between two persons
  p_trans <- 0.01                # Transmission proba when there is contact
  p_death <- 0.05                # Death proba

  max_time <- 5                  # Time for stopping the chain

  ## Theoretical mean number of infected
  m <- (1 - p_death) * (n_pop * p_contact * p_trans + 1)

  ## As functions for simus
  p_max_fct <- function(x){rep(p_trans,length(x))}
  p_Exit_fct  <- function(x){rep(p_death,length(x))}
  proba <- function(t, p_max){p_max}
  time_contact = function(x){rbinom(x, n_pop, p_contact)}

  ## Simulations
  nrep <- 1000
  size_pop <- matrix(NA_integer_, nrep, max_time)

  set.seed(20110721)             # A tout seigneur tout honneur
  for (rep in 1:nrep) {
    # Sim
    test.nosoi <- suppressMessages(
      nosoiSim(type="single", popStructure="none",
               length = max_time,
               max.infected = 10000,
               init.individuals = 1,
               nContact = time_contact,
               pTrans = proba,
               param.pTrans = list(p_max=p_max_fct),
               pExit = p_Exit_fct,
               param.pExit = NA
      )
    )
    # Save active pop
    size_pop[rep, ] <- sapply(1:max_time,
                              function(x) numberInfectedBGW(test.nosoi, x))
    # Print progress
    if (rep %% 100 == 0) print(rep)
  }

  empirical_mean_size <- colMeans(size_pop)

  theo_exp <- m ^ (1:max_time)

  expect_equal(theo_exp, empirical_mean_size, tolerance = 0.1)

})

test_that("Distribution of the number of active individuals - subcritic", {
  ## Parameters
  n_pop <- 100000                # Population size
  p_contact <- 0.001             # Contact proba between two persons
  p_trans <- 0.001               # Transmission proba when there is contact
  p_death <- 0.3                 # Death proba

  max_time <- 5                  # Time for stopping the chain

  ## Theoretical mean number of infected
  m <- (1 - p_death) * (n_pop * p_contact * p_trans + 1)

  ## As functions for simus
  p_max_fct <- function(x){rep(p_trans,length(x))}
  p_Exit_fct  <- function(x){rep(p_death,length(x))}
  proba <- function(t, p_max){p_max}
  time_contact = function(x){rbinom(x, n_pop, p_contact)}

  ## Simulations
  nrep <- 1000
  size_pop <- matrix(NA_integer_, nrep, max_time)

  set.seed(20110721)             # A tout seigneur tout honneur
  for (rep in 1:nrep) {
    # Sim
    test.nosoi <- suppressMessages(
      nosoiSim(type="single", popStructure="none",
               length = max_time,
               max.infected = 10000,
               init.individuals = 1,
               nContact = time_contact,
               pTrans = proba,
               param.pTrans = list(p_max=p_max_fct),
               pExit = p_Exit_fct,
               param.pExit = NA
      )
    )
    # Save active pop
    size_pop[rep, ] <- sapply(1:max_time,
                              function(x) numberInfectedBGW(test.nosoi, x))
    # Print progress
    if (rep %% 100 == 0) print(rep)
  }

  empirical_mean_size <- colMeans(size_pop)

  theo_exp <- m ^ (1:max_time)

  expect_equal(theo_exp, empirical_mean_size, tolerance = 0.1)

})
