#################################################################################################
#                                      IMPORT PACKAGES                                          #
#################################################################################################

#library(plyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(grid)
library(raster)
library(foreach)
library(doParallel)
library(readr)
library(dplyr)

setwd("~/Dropbox/C-Leuven Projects/ArboSim/ArboSim GeoEpi/discrete")

#################################################################################################
#                                          FUNCTIONS                                            #
#################################################################################################

#My formula:
viral.load.mos = function(x, Pmax, EIP, slope) {
  Pmax*1/2*(1+tanh((x-EIP)/slope))
}

#Albin's formula:
viral.load.mos = function(x, Pmax, EIP, slope) {
  Pmax*1/2*(1+tanh((x-EIP)/slope))}

viremia = function(x, Pmax, IIP, slopeIIP, Recov, slopeRecov) {
  Pmax*1/2*((1+tanh((x-IIP)/slopeIIP))-(1+tanh((x-Recov)/slopeRecov)))
}

viremia2 = function(x, Pmax, IIP, slopeIIP, Recov, slopeRecov, dist) {
  (Pmax*1/2*((1+tanh((x-IIP)/slopeIIP))-(1+tanh((x-Recov)/slopeRecov))))-dist
}

rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

#################################################################################################
#                                     MODEL PARAMETERS                                          #
#################################################################################################

# Transition matrix ---------------------------------------------------------------------------

AirlineMatrix <- read_delim("data/AirlineMatrix.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
AirlineData = melt(AirlineMatrix, id.variable="X1")
colnames(AirlineData) = c("From", "To", "Seats")

countries=unique(AirlineData$From)

AirlineData = merge(AirlineData,(AirlineData %>% group_by(From) %>% summarise(N=sum(Seats)))) %>% mutate(prob=Seats/N)

countries_centroid <- read_delim("data/countries_centroid.csv", 
                                 ";", escape_double = FALSE, trim_ws = TRUE)

countries_centroid = subset(countries_centroid, Code2 %in% countries)

countries_centroid = subset(countries_centroid, Country != "USA")

#Mosquito suitability 

AedesSuitability <- read_delim("data/AedesSuitability.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
#Suitability index will influence probability with which a human will get biten daily by an Aedes mosquito
correc.factor = 10 #because value are quite low, i'll multiply by 10

human.daily.biting.prob = AedesSuitability %>% mutate(prob=Aedes*10)

# Daily travel probability --------------------------------------------------------------------

Daily.travel.prob = 1/100

individuals=1 #number of initially infected individuals

#Dynamics of viremia in vertebrate (see EIP paper for justification)
IIPmean=4.844
IIPsd=0.5
slopeIIP=2.536
slopeIIP.sd=0.2
Pmax.vert.putative=1
Pmax.vert.var=0.002
Pmax.vert=Pmax.vert.putative-Pmax.vert.var
#Beta law parameters
alpha.vert=(Pmax.vert*(1-Pmax.vert)/Pmax.vert.var)*Pmax.vert
beta.vert=(alpha.vert/Pmax.vert)-alpha.vert
Recov=8.392
Recov.sd=0.5
slopeRecov=4.038
slopeRecov.sd=0.2
Prob.asympto = 2/3

individuals=1 #initially infected

#Dynamics of viremia in mosquitoes
EIPmean <- 7
slopeMos <-2.5
Pmax.mos.putative <- 0.8


surv.rate.vector=0.9 #probability of a vector to die each day
biting.rate.vector=0.3 #daily biting rate of a vector

#################################################################################################
#                                 PREPARE INITIAL HOSTS                                         #
#################################################################################################

hosts.ID =paste("H",c(1:individuals),sep="-")
inf.by = as.character(rep(NA,individuals))
attempt <-10
for(i in 1:attempt){
  table1.vert = data.table(hosts.ID,inf.by)
  table1.vert$IIP50 = rnorm(individuals,IIPmean,IIPsd)
  table1.vert$slopeInfec = rnorm(individuals,slopeIIP,slopeIIP.sd)
  table1.vert$Pmax = rbeta(individuals,alpha.vert,beta.vert)
  table1.vert$slopeRecov = rnorm(individuals,slopeRecov,slopeRecov.sd)
  table1.vert$Infec.pop = as.numeric(rep(NA,individuals))
  table1.vert$inf.in = sample(countries,individuals)
  table1.vert$inf.date = as.numeric(rep(0,individuals))
  table1.vert$asymptomatic = sample(c(TRUE, FALSE), individuals, replace=TRUE, prob=c(Prob.asympto, 1-Prob.asympto))
  table1.vert$active = TRUE
  
  table1.vert = table1.vert %>% group_by(hosts.ID) %>% mutate(Recov50=rtnorm(1,Recov,Recov.sd, a=IIP50+1, b=Inf))
  
  table1.vert = tryCatch(table1.vert %>% group_by(hosts.ID) %>% mutate(End.date=uniroot(viremia2, Pmax=Pmax, IIP=IIP50, slopeIIP=slopeInfec, Recov=Recov50, slopeRecov=slopeRecov, dist=0.01, upper=Recov50*2, lower=Recov50, extendInt="yes", maxiter=100000)$root), error=function(e) "skip")
  if(is.data.frame(table1.vert)==TRUE){break}
}
table1.vert = as.data.table(table1.vert)
setkey(table1.vert, hosts.ID)

H.count=individuals

####Movement history####
table1.mov = data.table(hosts.ID,origin=table1.vert[1:individuals,inf.in],now=table1.vert[1:individuals,inf.in],t=0)
setkey(table1.mov, hosts.ID)

mov.archiv=table1.mov

#################################################################################################
#                                 PREPARE MOSQUITOES HOSTS                                      #
#################################################################################################
#Switch version
table1.vect = as.data.table(NULL)
active.mosquitoes = NULL
M.count=0
#################################################################################################
#                                  SIMULATION (tests)                                           #
#################################################################################################

sim.length = 400

for (t in 1:sim.length){
  
#Step 1: Active hosts & dying mosquitoes ------------------------------------------------------------------------
  
  if(length(table1.vect) > 0){
  #Mosquitoes dying
  surviving = sample(c(TRUE,FALSE),length(active.mosquitoes),replace=TRUE,prob=c(surv.rate.vector,1-surv.rate.vector))
  table1.vect[as.character(active.mosquitoes[!surviving]),dead.date:=as.numeric(t)]
  table1.vect[as.character(active.mosquitoes[!surviving]),alive:=FALSE]

  active.mosquitoes = subset(table1.vect, alive==TRUE)$ID.vect #active mosquitoes
  }
  
  #Active hosts at this time point
  
  active.hosts = subset(subset(table1.vert,active==TRUE), End.date+inf.date > t)$hosts.ID #active hosts  
  table1.vert[!as.character(active.hosts),active:=FALSE]

#if(length(active.hosts) == 0 & length(active.mosquitoes) == 0){break}


#Step 1: Active hosts may be moving ------------------------------------------------------------------------

  table1.mov=subset(table1.mov,hosts.ID %in% active.hosts)
  
    for (j in active.hosts){
      setkey(table1.mov,hosts.ID)
      table1.mov[j,origin:=table1.mov[j,now]]
      moving=sample(c(TRUE,FALSE), 1, prob=c(Daily.travel.prob,1-Daily.travel.prob))#is individual moving?
      
      if(moving==TRUE){ 
        table1.mov[j,origin]
        
        going.to = sample(subset(AirlineData, From==table1.mov[j,origin])$To,1,replace=TRUE,prob=subset(AirlineData, From==table1.mov[j,origin])$prob)
        
        table1.mov[j,now:=going.to]
      }else{ 
        staying = table1.mov[j,origin]
        table1.mov[j,now:=staying]
      }
    } 

table1.mov$t = t
mov.archiv = plyr::rbind.fill(mov.archiv,table1.mov)

#Step 2: Active hosts may be biten in their new location ------------------------------------------------------------------------

for (j in active.hosts){
  bitten=sample(c(TRUE,FALSE), 1, prob=c(subset(human.daily.biting.prob, Location==table1.mov[j,now])$prob,1-subset(human.daily.biting.prob, Location==table1.mov[j,now])$prob))#is mosquito biting
  
  if(bitten == TRUE){
    #is bite infectious?
    x = t-table1.vert[j, inf.date]
    a.Pmax = table1.vert[j, Pmax]
    a.IIP50 =table1.vert[j, IIP50]
    a.slopeInfec =table1.vert[j, slopeInfec]
    a.Recov50 =table1.vert[j, Recov50]
    a.slopeRecov =table1.vert[j, slopeRecov]
    
    Ptransmit.t.vert2 = viremia(x,
                                a.Pmax,
                                a.IIP50,
                                a.slopeInfec ,
                                a.Recov50,
                                a.slopeRecov)
    
    Ptransmit.t.vert = ifelse(Ptransmit.t.vert2 > 0, Ptransmit.t.vert2, 0)
    transmitVert2Vect=sample(c(TRUE,FALSE),1,prob=c(Ptransmit.t.vert, 1-Ptransmit.t.vert))
    
    if(transmitVert2Vect==TRUE){ #mosquito gets infected
      
      ID.vect = as.character(paste("M",table1.mov[j,now],M.count+1,sep="-"))
      M.count = M.count+1
      
      table1.ongoing = data.frame(ID.vect)
      table1.ongoing$loc = table1.mov[j,now]
      table1.ongoing$dead.date=as.numeric(NA)
      table1.ongoing$inf.by=j
      table1.ongoing$infected.date=as.numeric(t)
      table1.ongoing$Infec.pop = as.numeric(Ptransmit.t.vert)
      table1.ongoing$alive = TRUE
      
        u = rbinom(1,1, Pmax.mos.putative)
        if(u == 1){
          result = rlogis(1, location = EIPmean, scale = 1/slopeMos)
        }else{result = Inf}

      result[result < 0] = 0
      
      table1.ongoing$EIP50 = result
      table1.ongoing = as.data.table(table1.ongoing)
      table1.vect = rbind(table1.vect,table1.ongoing)
      setkey(table1.vect,ID.vect)
    }
  } 
  t
} 

#Step 3: Active mosquitoes may bite in their location ------------------------------------------------------------------------

for (j in active.mosquitoes){
  bitting=sample(c(TRUE,FALSE), 1, prob=c(biting.rate.vector,1-biting.rate.vector))#is mosquito biting
  
  if(bitting == TRUE){#is bite infectious?
    
    x = t-table1.vect[j, infected.date]
    b.EIP50 = table1.vect[j, EIP50]
    
    Ptransmit.t.vect = ifelse(x >= b.EIP50, 1, 0)
    transmitVect2Vert=sample(c(TRUE,FALSE),1,prob=c(Ptransmit.t.vect, 1-Ptransmit.t.vect))
    
    if(transmitVect2Vert == TRUE){#New host is infected
      
      #new host biology
      
      ID.ongoing = paste("H",H.count+1,sep="-")
      H.count = H.count+1
      table1.vert.ongoing = data.frame(hosts.ID=ID.ongoing)
      table1.vert.ongoing$inf.by = table1.vect[j, ID.vect]

      attempt <-10
      for(i in 1:attempt){
        table1.vert.ongoing$IIP50 = rnorm(1,IIPmean,IIPsd)
        table1.vert.ongoing$slopeInfec = rnorm(1,slopeIIP,slopeIIP.sd)
        table1.vert.ongoing$Pmax = rbeta(1,alpha.vert,beta.vert)
        table1.vert.ongoing$slopeRecov = rnorm(1,slopeRecov,slopeRecov.sd)
        table1.vert.ongoing$Infec.pop = as.numeric(rep(NA,1))
        table1.vert.ongoing$inf.in = table1.vect[j, loc]
        table1.vert.ongoing$inf.date = t
        table1.vert.ongoing$asymptomatic = sample(c(TRUE, FALSE), 1, replace=TRUE, prob=c(Prob.asympto, 1-Prob.asympto))
        table1.vert.ongoing$active = TRUE
          
        table1.vert.ongoing = table1.vert.ongoing %>% group_by(hosts.ID) %>% mutate(Recov50=rtnorm(1,Recov,Recov.sd, a=IIP50+1, b=Inf))
        
        table1.vert.ongoing = tryCatch(table1.vert.ongoing %>% group_by(hosts.ID) %>% mutate(End.date=uniroot(viremia2, Pmax=Pmax, IIP=IIP50, slopeIIP=slopeInfec, Recov=Recov50, slopeRecov=slopeRecov, dist=0.01, upper=Recov50*2, lower=Recov50, extendInt="yes", maxiter=100000)$root), error=function(e) "skip")
        if(is.data.frame(table1.vert)==TRUE){break}
      }
      table1.vert.ongoing = as.data.table(table1.vert.ongoing)

      
      table1.vert= rbind(table1.vert,table1.vert.ongoing)
      setkey(table1.vert, hosts.ID)
      #new host movement
      
      ####Movement history####
      table1.mov.ongoing = data.table(hosts.ID=ID.ongoing,origin=table1.vect[j, loc],now=table1.vect[j, loc],t=t)
      table1.mov=rbind(table1.mov,table1.mov.ongoing)
    }
  }
}
}
nrow(table1.vert)

############
write_csv(table1.vect,"vectors.csv")
write_csv(table1.vert,"humans.csv")
write_csv(mov.archiv,"humans-moves.csv")

vectors <- read_csv("vectors.csv")
humans <- read_csv("humans.csv")
mov <- read_csv("humans-moves.csv")

############ Data

countries_centroid
worldmap <- borders("world", colour="gray50", fill="#efede1",t=1:182) 

sim.length = 238

data.summary = NULL
table1.vect[is.na(table1.vect$dead.date),]$dead.date = 239


p1 <- ggplot(data=data.summary)
p1 <- p1 + geom_line(aes(x=t, y=(vert.infected.cumulative), color="Vertebrates"))
p1 <- p1 + geom_line(aes(x=t, y=(vect.infected.cumulative), color="Vectors"))
p1 <- p1 + theme_bw() + labs(x="Time (in days)",y="Cumulative number of infected individuals")

p2 <- ggplot(data=data.summary)
p2 <- p2 + geom_line(aes(x=t, y=(vert.infected.active), color="Vertebrates"))
p2 <- p2 + geom_line(aes(x=t, y=(vect.infected.active), color="Vectors"))
p2 <- p2 + theme_bw() + labs(x="Time (in days)",y="Number of infectious individuals")

######## MAP
 data.airlines2 = merge(AirlineData, countries_centroid, by.x="From", by.y="Code2")
 data.airlines3 = merge(data.airlines2, countries_centroid, by.x="To", by.y="Code2", suffixes=c("FR","TO"))

data.airlines3 = subset(data.airlines3, Seats > 0)

 library(ggplot2)
 library(gganimate)
library(viridis)

summed.mov = mov %>% group_by(t,now) %>% summarise(N=length(hosts.ID))
data2 = merge(summed.mov, countries_centroid, by.x="now", by.y="Code2")

ggplot() + 
  theme_dark() + worldmap2 +
geom_point(data=data2,aes(x=Long,y=Lat,size=as.numeric(N)), show.legend = FALSE, color="orange1") +
  scale_size_continuous(c(5,20)) + 
  #  geom_curve(data=data.airlines3, aes(x = LongFR, y = LatFR,
  #                                      xend = LongTO, yend = LatTO,
  #                                      alpha=Seats), curvature = -0.2) +
  coord_cartesian(xlim = c(min(data2$Long)-0.25,
                           -50),
                  ylim = c(min(data2$Lat)-0.25,
                           max(data2$Lat)+0.25)) +

  # Here comes the gganimate specific bits
 labs(title = "Time: {current_frame}") +
  transition_manual(t)

worldmap2 <- borders("world", colour="gray50", fill="black") 

ggplot() + 
  theme_dark() + worldmap2 +
   geom_curve(data=data.airlines3, aes(x = LongFR, y = LatFR,
                                        xend = LongTO, yend = LatTO,
                                        color=log10(Seats)), curvature = -0.2,alpha=0.4) +
  scale_color_viridis(option="magma") + 
  coord_cartesian(xlim = c(min(data2$Long)-0.25,
                           -50),
                  ylim = c(min(data2$Lat)-0.25,
                           max(data2$Lat)+0.25))

AedesSuitability2 = merge(AedesSuitability,countries_centroid,by.x="Location",by.y="Code2")

library(maps)
world <- map_data("world")
head(world)

AedesSuitability2[19,5] = "USA"

us_states_elec <- left_join(world, AedesSuitability2,by=c("region" = "Country"))

ggplot(data = us_states_elec,
            mapping = aes(x = long, y = lat,
                          group = group,order=order,fill=Aedes)) + geom_polygon() +
  scale_fill_viridis(option="magma") + 
  coord_cartesian(xlim = c(min(data2$Long)-0.25,
                           -50),
                  ylim = c(min(data2$Lat)-0.25,
                           max(data2$Lat)+0.25))
length(unique(us_states_elec$Aedes))
