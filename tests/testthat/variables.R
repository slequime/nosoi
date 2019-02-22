t_incub <- function(x){rnorm(x,mean = 5,sd=1)}
p_max <- function(x){rbeta(x,shape1 = 5,shape2=2)}

proba <- function(t,p_max,t_incub){
  if(t <= t_incub){p=0}
  if(t >= t_incub){p=p_max}
  return(p)
}

time_contact = function(x){round(rnorm(x, 3, 1), 0)}
