## Earth has surpassed its sustainable human carrying capacity
## Corey Bradshaw
## Flinders University, Adelaide, South Australia
## Global Ecology Laboratory | Partuyarta Ngadluku Wardli Kuu
## February 2024 / updated February 2026

# required R libraries
library(boot)
library(dismo)
library(gbm)
library(ggplot2)
library(gridExtra)
library(lmtest)
library(performance)
library(plotrix)
library(sjPlot)
library(sandwich)
library(tmvnsim)
library(truncnorm)
library(wCorr)

# source files
source("new_lmer_AIC_tables3.R")
source("r.squared.R")

# functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

# geometric mean
gmMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## sample avoiding consecutive values of order o
NonSeqSample <- function(x, size, ord=1, replace) {
  # x = vector to sample
  # size = size of resample
  # ord = order to consider (1 = no consecutive, 2 = no 2nd-order consecutive ...)
  # replace = with (TRUE) or without (FALSE) replacement
  vsmp <- sort(sample(x=x, size=size, replace=F))
  diff.vsmp <- diff(vsmp)
  vsmp <- vsmp[-(which(diff.vsmp <= ord)+1)]
  new.size <- size - length(vsmp)
  
  while(new.size >= 1) {
    vsmp.add <- sort(sample(x=x, size=new.size+1, replace=F))
    vsmp <- sort(c(vsmp, vsmp.add))
    diff.vsmp <- diff(vsmp)
    vsmp <- vsmp[-(which(diff.vsmp <= ord)+1)]
    new.size <- size - length(vsmp)
  }
  return(vsmp)
}

## Historical estimates of world population
## census.gov/data/tables/time-series/demo/international-programs/historical-est-worldpop.html
hpop <- read.csv('Npre1950.csv', sep=",", header=T)
head(hpop)
hpop$pop.md <- apply(hpop[,c(2,3)], MARGIN=1, median, na.rm=T)
head(hpop)
hpop.out <- data.frame(hpop[,1], hpop[,4], hpop[,3], hpop[,2])
colnames(hpop.out) <- c("year", "popMD", "popUP", "popLO")
head(hpop.out)

## UN population
## https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1
UNpop <- read.csv('UNpop.csv', sep=",", header=T)
head(UNpop)
UNtot <- subset(UNpop, Area=="Total" & Sex=="Both Sexes")
head(UNtot)
UNtot2022 <- subset(UNtot, Year==2022)
head(UNtot2022)
sum(UNtot2022$Value)/10^9

popdat <- read.csv("worldpophist.csv")
head(popdat)
tail(popdat)
popdat$r <- c(NA, log(popdat$pop[2:length(popdat$pop)] / popdat$pop[1:(length(popdat$pop)-1)]))
head(popdat)
popdat$rpcap <- popdat$r/popdat$pop
head(popdat)

par(mfrow=c(1,3))
plot(popdat$year, popdat$pop, type="l")
plot(popdat$year, popdat$r, type="l")
plot(popdat$year, popdat$rpcap, type="l")
par(mfrow=c(1,1))

## export popdat
write.csv(popdat, "pop_r_out.csv", row.names=F)

## interpolate yearly
yrintp <- seq(-10000, 2023, 1)
popintp <- approx(popdat$year, popdat$pop, xout = yrintp)

popdatintp <- data.frame(yrintp, popintp$y)
colnames(popdatintp) <- c("year", "pop")
popdatintp$r <- c(NA, log(popdatintp$pop[2:length(popdatintp$pop)] / popdatintp$pop[1:(length(popdatintp$pop)-1)]))
popdatintp$rpcap <- popdatintp$r/popdatintp$pop
head(popdatintp)

popdatintp$popdiff <- c(NA,diff(popdatintp$pop))
plot(popdatintp$year[2:(length(popdatintp$year))], (diff(popdatintp$pop)), type="l")

## assume different errors on population estimates
head(popdat)

## split into 3 phases
popPhase1 <- subset(popdat, year >= 1800 & year < 1950)
popPhase2 <- subset(popdat, year >= 1950 & year < 1962)
popPhase3 <- subset(popdat, year >= 1962)

## set uncertainties
pcUncertPhase1 <- 0.05 # 5% uncertainty
pcUncertPhase2 <- 0.02 # 2% uncertainty
pcUncertPhase3 <- 0.01 # 1% uncertainty


## Ricker and Gompertz models for each phase and examine relative AIC weights,
## model-averaged K estimates, and confidence intervals on r-vs-Nt slopes

## Phase 3
## Ricker
plot((popPhase3$pop[1:(length(popPhase3$pop)-1)]), popPhase3$r[2:length(popPhase3$pop)], pch=19, xlab="Nt", ylab="r")
y.R.Phase3 <- popPhase3$r[2:length(popPhase3$pop)]
x.R.Phase3 <- popPhase3$pop[1:(length(popPhase3$pop)-1)]
dat.R.Phase3 <- data.frame(x.R.Phase3, y.R.Phase3)
head(dat.R.Phase3)
fitlin3 <- lm(y.R.Phase3 ~ x.R.Phase3, data=dat.R.Phase3)
summary(fitlin3)
abline(fitlin3, lty=2, col="red")

## calculate 95% prediction confidence interval
x.pred.R.Phase3 <- data.frame(x.R.Phase3 = seq(min(dat.R.Phase3$x.R.Phase3), max(dat.R.Phase3$x.R.Phase3), length.out=100))
x.pred.R.Phase3 <- data.frame(x.R.Phase3 = seq(0, 20e9, length.out=100))
pred.interval.R.Phase3 <- predict(fitlin3, newdata=x.pred.R.Phase3, interval="confidence", level=0.95)
lines(x.pred.R.Phase3$x.R.Phase3, pred.interval.R.Phase3[,"lwr"], col="blue", lty=2)
lines(x.pred.R.Phase3$x.R.Phase3, pred.interval.R.Phase3[,"upr"], col="blue", lty=2)

k3 <- as.numeric(-coef(fitlin3)[1]/coef(fitlin3)[2])
round(k3/10^9, 4)

# confidence interval for k
fitlin3.summ <- summary(fitlin3)
fitlin3.summ$coefficients

int3.up <- confint(fitlin3, level=0.95)[1,2]
int3.lo <- confint(fitlin3, level=0.95)[1,1]
print(c(int3.lo, int3.up))

slp3.up <- confint(fitlin3, level=0.95)[2,2]
slp3.lo <- confint(fitlin3, level=0.95)[2,1]
print(c(slp3.lo, slp3.up))

k3.lo <- as.numeric(-int3.up/slp3.lo)
k3.up <- as.numeric(-int3.lo/slp3.up)
print(c(k3.lo, k3, k3.up)/10^9)

# calculate time to stability (where r = 0)
pop.start <- popPhase3$pop[which(popPhase3$year==1962)]
pop.end <- popPhase3$pop[which(popPhase3$year==2023)]
r.pred.start <- predict(fitlin3, newdata=data.frame(x.R.Phase3=popPhase3$pop[which(popPhase3$year==1962)]))
r.pred.end <- predict(fitlin3, newdata=data.frame(x.R.Phase3=popPhase3$pop[which(popPhase3$year==2023)]))
year.start <- popPhase3$year[1]
year.end <- popPhase3$year[length(popPhase3$year)]

pastC2 <- sqrt((pop.end - pop.start)^2 + (r.pred.start - r.pred.end)^2)
pastC2rate <- pastC2 / (year.end-year.start)

futH.lo <- sqrt((k3.lo - pop.end)^2 + (r.pred.end - 0)^2)
stabilityT.lo <- round(futH.lo/pastC2rate, 0)
year.end + stabilityT.lo

futH.up <- sqrt((k3.up - pop.end)^2 + (r.pred.end - 0)^2)
stabilityT.up <- round(futH.up/pastC2rate, 0)
year.end + stabilityT.up

## get intercept range (r)
int3.R <- as.numeric(coef(fitlin3)[1])
int3.R.lo <- confint(fitlin3, level=0.95)[1,1]
int3.R.up <- confint(fitlin3, level=0.95)[1,2]
print(c(int3.R.lo, int3.R, int3.R.up))

## Gompertz
plot(log(popPhase3$pop[1:(length(popPhase3$pop)-1)]), popPhase3$r[2:length(popPhase3$pop)], pch=19, xlab="log Nt", ylab="r")
y.G.Phase3 <- popPhase3$r[2:length(popPhase3$pop)]
x.G.Phase3 <- log(popPhase3$pop[1:(length(popPhase3$pop)-1)])
dat.G.Phase3 <- data.frame(x.G.Phase3, y.G.Phase3)
head(dat.G.Phase3)
fitlog3 <- lm(y.G.Phase3 ~ x.G.Phase3, data=dat.G.Phase3)
abline(fitlog3, lty=2, col="red")

## calculate 95% prediction confidence interval
x.pred.G.Phase3 <- data.frame(x.G.Phase3 = seq(min(dat.G.Phase3$x.G.Phase3), max(dat.G.Phase3$x.G.Phase3), length.out=100))
x.pred.G.Phase3 <- data.frame(x.G.Phase3 = log(seq(0, 20e9, length.out=100)))
pred.interval.G.Phase3 <- predict(fitlog3, newdata=x.pred.G.Phase3, interval="confidence", level=0.95)
lines(x.pred.G.Phase3$x.G.Phase3, pred.interval.G.Phase3[,"lwr"], col="blue", lty=2)
lines(x.pred.G.Phase3$x.G.Phase3, pred.interval.G.Phase3[,"upr"], col="blue", lty=2)

klog3 <- exp(as.numeric(-coef(fitlog3)[1]/coef(fitlog3)[2]))
round(klog3/10^9, 4)

Phase3.results.out <- data.frame(x.R=x.pred.R.Phase3/10^9,
                                 R.fit=pred.interval.R.Phase3[,"fit"],
                                 R.upr=pred.interval.R.Phase3[,"upr"],
                                 R.lwr=pred.interval.R.Phase3[,"lwr"],
                                 x.G=exp(x.pred.G.Phase3)/10^9,
                                 G.fit=pred.interval.G.Phase3[,"fit"],
                                 G.upr=pred.interval.G.Phase3[,"upr"],
                                 G.lwr=pred.interval.G.Phase3[,"lwr"])
head(Phase3.results.out)

# save to .csv
write.csv(Phase3.results.out, "Phase3_Ricker_Gompertz_fits.csv", row.names=F)

## plot both Ricker and Gompertz fits on linear scale in ggplot
ggplot() +
  geom_point(data=dat.R.Phase3, aes(x=x.R.Phase3, y=y.R.Phase3), color="black") +
  geom_line(aes(x=x.pred.R.Phase3$x.R.Phase3, y=pred.interval.R.Phase3[,"fit"]), color="red", linetype="dashed") +
  geom_ribbon(aes(x=x.pred.R.Phase3$x.R.Phase3, ymin=pred.interval.R.Phase3[,"lwr"], ymax=pred.interval.R.Phase3[,"upr"]), alpha=0.2, fill="blue") +
  geom_line(aes(x=exp(x.pred.G.Phase3$x.G.Phase3), y=pred.interval.G.Phase3[,"fit"]), color="darkgreen", linetype="dashed") +
  geom_ribbon(aes(x=exp(x.pred.G.Phase3$x.G.Phase3), ymin=pred.interval.G.Phase3[,"lwr"], ymax=pred.interval.G.Phase3[,"upr"]), alpha=0.2, fill="lightgreen") +
  labs(x="Nt", y="r") +
  xlim(min(dat.R.Phase3$x.R.Phase3), max(dat.R.Phase3$x.R.Phase3)) +
  ylim(min(dat.R.Phase3$y.R.Phase3), max(dat.R.Phase3$y.R.Phase3)) +
  theme_minimal()

RLGL3.AIC.vec <- c(AICc(fitlin3), AICc(fitlog3))
RLGL3.dAIC.vec <- delta.IC(RLGL3.AIC.vec)
RLGL3.wAIC.vec <- weight.IC(RLGL3.dAIC.vec)
print(RLGL3.wAIC.vec)

ER.RLGL3 <- RLGL3.wAIC.vec[1]/RLGL3.wAIC.vec[2]
ER.RLGL3

model.avg.Kmax3 <- round(sum(RLGL3.wAIC.vec * c(k3, klog3)) / 10^9, 4)
model.avg.Kmax3


## Phase 1
## Ricker
y.R.Phase1 <- popPhase1$r[2:length(popPhase1$pop)]
x.R.Phase1 <- popPhase1$pop[1:(length(popPhase1$pop)-1)]
dat.R.Phase1 <- data.frame(x.R.Phase1, y.R.Phase1)
head(dat.R.Phase1)

plot(x.R.Phase1, y.R.Phase1, pch=19, xlab="Nt", ylab="r")
fitlin1 <- lm(y.R.Phase1 ~ x.R.Phase1)
abline(fitlin1, lty=2, col="red")
summary(fitlin1)

## calculate 95% prediction confidence interval
x.pred.R.Phase1 <- data.frame(x.R.Phase1 = seq(min(dat.R.Phase1$x.R.Phase1), max(dat.R.Phase1$x.R.Phase1), length.out=100))
x.pred.R.Phase1 <- data.frame(x.R.Phase1 = seq(0, 2.47464809e9, length.out=100))
pred.interval.R.Phase1 <- predict(fitlin1, newdata=x.pred.R.Phase1, interval="confidence", level=0.95)
lines(x.pred.R.Phase1$x.R.Phase1, pred.interval.R.Phase1[,"lwr"], col="blue", lty=2)
lines(x.pred.R.Phase1$x.R.Phase1, pred.interval.R.Phase1[,"upr"], col="blue", lty=2)

## Gompertz
y.G.Phase1 <- popPhase1$r[2:length(popPhase1$pop)]
x.G.Phase1 <- log(popPhase1$pop[1:(length(popPhase1$pop)-1)])
dat.G.Phase1 <- data.frame(x.G.Phase1, y.G.Phase1)
head(dat.G.Phase1)

plot(x.G.Phase1, y.G.Phase1, pch=19, xlab="log Nt", ylab="r")
fitlog1 <- lm(y.G.Phase1 ~ x.G.Phase1)
abline(fitlog1, lty=2, col="red")

## calculate 95% prediction confidence interval
x.pred.G.Phase1 <- data.frame(x.G.Phase1 = seq(min(dat.G.Phase1$x.G.Phase1), max(dat.G.Phase1$x.G.Phase1), length.out=100))
x.pred.G.Phase1 <- data.frame(x.G.Phase1 = log(seq(0, 2.47464809e9, length.out=100)))
pred.interval.G.Phase1 <- predict(fitlog1, newdata=x.pred.G.Phase1, interval="confidence", level=0.95)
lines(x.pred.G.Phase1$x.G.Phase1, pred.interval.G.Phase1[,"lwr"], col="blue", lty=2)
lines(x.pred.G.Phase1$x.G.Phase1, pred.interval.G.Phase1[,"upr"], col="blue", lty=2)

# confidence interval for k
fitlin1.summ <- summary(fitlin1)
fitlin1.summ$coefficients

int1.up <- confint(fitlin1, level=0.95)[1,2]
int1.lo <- confint(fitlin1, level=0.95)[1,1]
print(c(int1.lo, int1.up))

slp1.up <- confint(fitlin1, level=0.95)[2,2]
slp1.lo <- confint(fitlin1, level=0.95)[2,1]
print(c(slp1.lo, slp1.up))


## get intercept range (r)
int1.R <- as.numeric(coef(fitlin1)[1])
int1.R.lo <- confint(fitlin1, level=0.95)[1,1]
int1.R.up <- confint(fitlin1, level=0.95)[1,2]
print(c(int1.R.lo, int1.R, int1.R.up))

int1.G <- as.numeric(coef(fitlog1)[1])
int1.G.lo <- confint(fitlog1, level=0.95)[1,1]
int1.G.up <- confint(fitlog1, level=0.95)[1,2]
print(c(int1.G.lo, int1.G, int1.G.up))

## export results
Phase1.results.out <- data.frame(x.R=x.pred.R.Phase1/10^9,
                                 R.fit=pred.interval.R.Phase1[,"fit"],
                                 R.upr=pred.interval.R.Phase1[,"upr"],
                                 R.lwr=pred.interval.R.Phase1[,"lwr"],
                                 x.G=exp(x.pred.G.Phase1)/10^9,
                                 G.fit=pred.interval.G.Phase1[,"fit"],
                                 G.upr=pred.interval.G.Phase1[,"upr"],
                                 G.lwr=pred.interval.G.Phase1[,"lwr"])
head(Phase1.results.out)

# save to .csv
write.csv(Phase1.results.out, "Phase1_Ricker_Gompertz_fits.csv", row.names=F)


RLGL1.AIC.vec <- c(AICc(fitlin1), AICc(fitlog1))
RLGL1.dAIC.vec <- delta.IC(RLGL1.AIC.vec)
RLGL1.wAIC.vec <- weight.IC(RLGL1.dAIC.vec)
print(RLGL1.wAIC.vec)
ER.RLGL1 <- RLGL1.wAIC.vec[1]/RLGL1.wAIC.vec[2]
ER.RLGL1

## K (Nmax) for Phase 1
K.phase1 <- max(dat.R.Phase1$x.R.Phase1) / 1e9
K.phase1

## plot both Ricker and Gompertz on linear scale in ggplot
ggplot() +
  geom_point(data=dat.R.Phase1, aes(x=x.R.Phase1, y=y.R.Phase1), color="black") +
  geom_line(aes(x=x.pred.R.Phase1$x.R.Phase1, y=pred.interval.R.Phase1[,"fit"]), color="red", linetype="dashed") +
  geom_ribbon(aes(x=x.pred.R.Phase1$x.R.Phase1, ymin=pred.interval.R.Phase1[,"lwr"], ymax=pred.interval.R.Phase1[,"upr"]), alpha=0.2, fill="blue") +
  geom_line(aes(x=exp(x.pred.G.Phase1$x.G.Phase1), y=pred.interval.G.Phase1[,"fit"]), color="darkgreen", linetype="dashed") +
  geom_ribbon(aes(x=exp(x.pred.G.Phase1$x.G.Phase1), ymin=pred.interval.G.Phase1[,"lwr"], ymax=pred.interval.G.Phase1[,"upr"]), alpha=0.2, fill="lightgreen") +
  labs(x="Nt", y="r") +
  xlim(min(dat.R.Phase1$x.R.Phase1), max(dat.R.Phase1$x.R.Phase1)) +
  ylim(min(dat.R.Phase1$y.R.Phase1), max(dat.R.Phase1$y.R.Phase1)) +
  theme_minimal()


## Phase 2
## Ricker
# remove 1949-1950 transition
popPhase2.clean <- popPhase2[-which(popPhase2$year == 1950), ]
popPhase2.clean

y.R.Phase2 <- popPhase2.clean$r[2:length(popPhase2.clean$pop)]
x.R.Phase2 <- popPhase2.clean$pop[1:(length(popPhase2.clean$pop)-1)]
dat.R.Phase2 <- data.frame(x.R.Phase2, y.R.Phase2)
head(dat.R.Phase2)

plot(x.R.Phase2, y.R.Phase2, pch=19, xlab="Nt", ylab="r")
fitlin2 <- lm(y.R.Phase2 ~ x.R.Phase2)
abline(fitlin2, lty=2, col="red")
summary(fitlin2)

## calculate 95% prediction confidence interval
x.pred.R.Phase2 <- data.frame(x.R.Phase2 = seq(min(dat.R.Phase2$x.R.Phase2), max(dat.R.Phase2$x.R.Phase2), length.out=100))
x.pred.R.Phase2 <- data.frame(x.R.Phase2 = seq(0, 3e9, length.out=100))
pred.interval.R.Phase2 <- predict(fitlin2, newdata=x.pred.R.Phase2, interval="confidence", level=0.95)
lines(x.pred.R.Phase2$x.R.Phase2, pred.interval.R.Phase2[,"lwr"], col="blue", lty=2)
lines(x.pred.R.Phase2$x.R.Phase2, pred.interval.R.Phase2[,"upr"], col="blue", lty=2)

## Gompertz
y.G.Phase2 <- popPhase2$r[2:length(popPhase2$pop)]
x.G.Phase2 <- log(popPhase2$pop[1:(length(popPhase2$pop)-1)])
dat.G.Phase2 <- data.frame(x.G.Phase2, y.G.Phase2)
head(dat.G.Phase2)

plot(x.G.Phase2, y.G.Phase2, pch=19, xlab="log Nt", ylab="r")
fitlog2 <- lm(y.G.Phase2 ~ x.G.Phase2)
abline(fitlog2, lty=2, col="red")

## calculate 95% prediction confidence interval
x.pred.G.Phase2 <- data.frame(x.G.Phase2 = seq(min(dat.G.Phase2$x.G.Phase2), max(dat.G.Phase2$x.G.Phase2), length.out=100))
x.pred.G.Phase2 <- data.frame(x.G.Phase2 = log(seq(0, 3e9, length.out=100)))
pred.interval.G.Phase2 <- predict(fitlog2, newdata=x.pred.G.Phase2, interval="confidence", level=0.95)
lines(x.pred.G.Phase2$x.G.Phase2, pred.interval.G.Phase2[,"lwr"], col="blue", lty=2)
lines(x.pred.G.Phase2$x.G.Phase2, pred.interval.G.Phase2[,"upr"], col="blue", lty=2)

## get intercept range (r)
int2.R <- as.numeric(coef(fitlin2)[1])
int2.R.lo <- confint(fitlin2, level=0.95)[1,1]
int2.R.up <- confint(fitlin2, level=0.95)[1,2]
print(c(int2.R.lo, int2.R, int2.R.up))

slope2.R <- as.numeric(coef(fitlin2)[2])
slope2.R.lo <- confint(fitlin2, level=0.95)[2,1]
slope2.R.up <- confint(fitlin2, level=0.95)[2,2]
print(c(slope2.R.lo, slope2.R, slope2.R.up))
print((c(slope2.R.lo*1e9, slope2.R*1e9, slope2.R.up*1e9)))

int2.G <- as.numeric(coef(fitlog2)[1])
int2.G.lo <- confint(fitlog2, level=0.95)[1,1]
int2.G.up <- confint(fitlog2, level=0.95)[1,2]
print(c(int2.G.lo, int2.G, int2.G.up))

## export results
Phase2.results.out <- data.frame(x.R=x.pred.R.Phase2/10^9,
                                 R.fit=pred.interval.R.Phase2[,"fit"],
                                 R.upr=pred.interval.R.Phase2[,"upr"],
                                 R.lwr=pred.interval.R.Phase2[,"lwr"],
                                 x.G=exp(x.pred.G.Phase2)/10^9,
                                 G.fit=pred.interval.G.Phase2[,"fit"],
                                 G.upr=pred.interval.G.Phase2[,"upr"],
                                 G.lwr=pred.interval.G.Phase2[,"lwr"])
head(Phase2.results.out)

# save to .csv
write.csv(Phase2.results.out, "Phase2_Ricker_Gompertz_fits.csv", row.names=F)


RLGL2.AIC.vec <- c(AICc(fitlin2), AICc(fitlog2))
RLGL2.dAIC.vec <- delta.IC(RLGL2.AIC.vec)
RLGL2.wAIC.vec <- weight.IC(RLGL2.dAIC.vec)
print(RLGL2.wAIC.vec)
ER.RLGL2 <- RLGL2.wAIC.vec[1]/RLGL2.wAIC.vec[2]
ER.RLGL2


## plot both Ricker and Gompertz on linear scale in ggplot
ggplot() +
  geom_point(data=dat.R.Phase2, aes(x=x.R.Phase2, y=y.R.Phase2), color="black") +
  geom_line(aes(x=x.pred.R.Phase2$x.R.Phase2, y=pred.interval.R.Phase2[,"fit"]), color="red", linetype="dashed") +
  geom_ribbon(aes(x=x.pred.R.Phase2$x.R.Phase2, ymin=pred.interval.R.Phase2[,"lwr"], ymax=pred.interval.R.Phase2[,"upr"]), alpha=0.2, fill="blue") +
  geom_line(aes(x=exp(x.pred.G.Phase2$x.G.Phase2), y=pred.interval.G.Phase2[,"fit"]), color="darkgreen", linetype="dashed") +
  geom_ribbon(aes(x=exp(x.pred.G.Phase2$x.G.Phase2), ymin=pred.interval.G.Phase2[,"lwr"], ymax=pred.interval.G.Phase2[,"upr"]), alpha=0.2, fill="lightgreen") +
  labs(x="Nt", y="r") +
  xlim(min(dat.R.Phase2$x.R.Phase2), max(dat.R.Phase2$x.R.Phase2)) +
  ylim(min(dat.R.Phase2$y.R.Phase2), max(dat.R.Phase2$y.R.Phase2)) +
  theme_minimal()


## autocorrelation in r series
plot_combined_acf_pacf <- function(ts_data, max_lag = 20) {
  
  # Calculate ACF and PACF
  acf_result <- acf(ts_data, lag.max = max_lag, plot = FALSE)
  pacf_result <- pacf(ts_data, lag.max = max_lag, plot = FALSE)
  
  # Confidence interval
  ci <- qnorm(0.975) / sqrt(length(ts_data))
  
  # Prepare ACF data
  acf_df <- data.frame(
    lag = acf_result$lag[-1],
    acf = acf_result$acf[-1]
  )
  
  # Prepare PACF data
  pacf_df <- data.frame(
    lag = pacf_result$lag,
    pacf = pacf_result$acf
  )
  
  # ACF plot
  p1 <- ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed", color = "blue") +
    geom_segment(aes(xend = lag, yend = 0)) +
    geom_point() +
    labs(x = "lag", y = "autocorrelation") +
    theme_minimal()
  
  # PACF plot
  p2 <- ggplot(pacf_df, aes(x = lag, y = pacf)) +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed", color = "blue") +
    geom_segment(aes(xend = lag, yend = 0)) +
    geom_point() +
    labs(x = "lag", y = "partial autocorrelation") +
    theme_minimal()
  
  # Combine plots
  grid.arrange(p1, p2, ncol = 1)
  
  # Return significant lags
  list(
    acf_significant = acf_df$lag[abs(acf_df$acf) > ci],
    pacf_significant = pacf_df$lag[abs(pacf_df$pacf) > ci]
  )
}

phase3.acf.pacf.results <- plot_combined_acf_pacf(popPhase3$r, max_lag = 40)
phase3.acf.pacf.results

phase1.acf.pacf.results <- plot_combined_acf_pacf(popPhase1$r, max_lag = 40)
phase1.acf.pacf.results

phase2.acf.pacf.results <- plot_combined_acf_pacf(popPhase2$r, max_lag = 40)
phase2.acf.pacf.results

## Durbin-Waatson test for autocorrelation
dwtest(fitlin3)

## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
# Phase 3
coeftest.Phase3 <- coeftest(fitlin3, vcov = NeweyWest(fitlin3, prewhite = FALSE))
coeftest.Phase3

coeft3.confint <- confint(coeftest.Phase3)
coeft3.confint
confint(fitlin3)

int3.NW.up <- coeft3.confint[1,2]
int3.NW.lo <- coeft3.confint[1,1]

slp3.NW.up <- coeft3.confint[2,2]
slp3.NW.lo <- coeft3.confint[2,1]
print(c(slp3.NW.lo, slp3.NW.up))*1e9

k3.NW.lo <- as.numeric(-int3.NW.up/slp3.NW.lo)
k3.NW.up <- as.numeric(-int3.NW.lo/slp3.NW.up)
print(c(k3.NW.lo, k3.NW.up)/10^9)


# Phase 1
coeftest.Phase1 <- coeftest(fitlin1, vcov = NeweyWest(fitlin1, prewhite = FALSE))
coeftest.Phase1

coeft1.confint <- confint(coeftest.Phase1)
coeft1.confint
confint(fitlin1)

int1.NW.up <- coeft1.confint[1,2]
int1.NW.lo <- coeft1.confint[1,1]

slp1.NW.up <- coeft1.confint[2,2]
slp1.NW.lo <- coeft1.confint[2,1]
print(c(slp1.NW.lo, slp1.NW.up))*1e9


# Phase 2
coeftest.Phase2 <- coeftest(fitlin2, vcov = NeweyWest(fitlin2, prewhite = FALSE))
coeftest.Phase2

coeft2.confint <- confint(coeftest.Phase2)
coeft2.confint
confint(fitlin2)

int2.NW.up <- coeft2.confint[1,2]
int2.NW.lo <- coeft2.confint[1,1]

slp2.NW.up <- coeft2.confint[2,2]
slp2.NW.lo <- coeft2.confint[2,1]
print(c(slp2.NW.lo, slp2.NW.up))*1e9


# calculate stochastic time series of r based on phase-specific certainties & estimate slope of
# relationship between r and Nt
iter <- 10000
Phase1rN.slope <- Phase2rN.slope <- Phase3rN.slope <- Phase3K <- rep(NA, iter)
for (i in 1:iter) {
  # Phase 1: facilitation
  Phase1pop.rsmp <- runif(length(popPhase1$pop), min = popPhase1$pop - (popPhase1$pop * pcUncertPhase1),
        max = popPhase1$pop + (popPhase1$pop * pcUncertPhase1))
  Phase1r.rsmp <- c(NA, log(Phase1pop.rsmp[2:length(Phase1pop.rsmp)] / Phase1pop.rsmp[1:(length(Phase1pop.rsmp)-1)]))
  Phase1rN.fit <- lm(Phase1r.rsmp ~ Phase1pop.rsmp)  
  #summary(Phase1rN.fit)
  Phase1rN.slope[i] <- as.numeric(coef(Phase1rN.fit)[2])
  
  # Phase 2: transition
  Phase2pop.rsmp <- runif(length(popPhase2$pop), min = popPhase2$pop - (popPhase2$pop * pcUncertPhase2),
                          max = popPhase2$pop + (popPhase2$pop * pcUncertPhase2))
  Phase2r.rsmp <- c(NA, log(Phase2pop.rsmp[2:length(Phase2pop.rsmp)] / Phase2pop.rsmp[1:(length(Phase2pop.rsmp)-1)]))
  Phase2rN.fit <- lm(Phase2r.rsmp ~ Phase2pop.rsmp)  
  #summary(Phase2rN.fit)
  Phase2rN.slope[i] <- as.numeric(coef(Phase2rN.fit)[2])
  
  # Phase 3: negative
  Phase3pop.rsmp <- runif(length(popPhase3$pop), min = popPhase3$pop - (popPhase3$pop * pcUncertPhase3),
                          max = popPhase3$pop + (popPhase3$pop * pcUncertPhase3))
  Phase3r.rsmp <- c(NA, log(Phase3pop.rsmp[2:length(Phase3pop.rsmp)] / Phase3pop.rsmp[1:(length(Phase3pop.rsmp)-1)]))
  Phase3rN.fit <- lm(Phase3r.rsmp ~ Phase3pop.rsmp)
  #summary(Phase3rN.fit)
  Phase3rN.slope[i] <- as.numeric(coef(Phase3rN.fit)[2])
  Phase3K[i] <- as.numeric(-coef(Phase3rN.fit)[1]/coef(Phase3rN.fit)[2]) / 10^9
  
}

# calculate 95% confidence intervals
Phase1rN.slope.lo <- quantile(Phase1rN.slope, probs=0.025, na.rm=T)
Phase1rN.slope.up <- quantile(Phase1rN.slope, probs=0.975, na.rm=T)
Phase2rN.slope.lo <- quantile(Phase2rN.slope, probs=0.025, na.rm=T)
Phase2rN.slope.up <- quantile(Phase2rN.slope, probs=0.975, na.rm=T)
Phase3rN.slope.lo <- quantile(Phase3rN.slope, probs=0.025, na.rm=T)
Phase3rN.slope.up <- quantile(Phase3rN.slope, probs=0.975, na.rm=T)

print(c(Phase1rN.slope.lo, Phase1rN.slope.up))*1e9
print(c(Phase2rN.slope.lo, Phase2rN.slope.up))*1e9
print(c(Phase3rN.slope.lo, Phase3rN.slope.up))*1e9

Phase3rN.K.lo <- quantile(Phase3K, probs=0.025, na.rm=T)
Phase3rN.K.up <- quantile(Phase3K, probs=0.975, na.rm=T)
print(c(Phase3rN.K.lo, Phase3rN.K.up))



####################################
## population x temperature anomaly 
## (metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.0.2.0/download.html)

popXta <- read.csv("popXtempanom.csv", header=T)
head(popXta)
tail(popXta)
popXta.pre1950 <- subset(popXta, year < 1950)
popXta.1950.1961 <- subset(popXta, year >= 1950 & year <= 1961)
popXta.post1961 <- subset(popXta, year > 1961)

#dat.use <- popXta.pre1950
#dat.use <- popXta.1950.1961
dat.use <- popXta.post1961[1:(length(popXta.post1961$year)-2), ]

iter <- 10000
itdiv <- iter/10

R2.vec <- ER.vec <- p.NW.vec <- slope.NW.vec <- R2.NW.vec <- rep(NA, iter)
for (i in 1:iter) {
  ta.samp <- runif(dim(dat.use)[1], dat.use$anomLO, dat.use$anomUP)
  R2.vec[i] <- linreg.ER(dat.use$pop, ta.samp)[2]
  ER.vec[i] <- linreg.ER(dat.use$pop, ta.samp)[1]
  
  popXta.fit <- lm(ta.samp ~ dat.use$pop)
  coeftest.fit <- coeftest(popXta.fit, vcov = NeweyWest(popXta.fit, prewhite = FALSE), save=T)
  
  coeftest.fitted <- attr(coeftest.fit, "object")[[5]]
  R2.NW.vec[i] <- cor(ta.samp, coeftest.fitted)^2
  slope.NW.vec[i] <- as.numeric(attr(coeftest.fit, "object")[[1]][2])
  p.NW.vec[i] <- coeftest.fit[7]
  
  if (i %% itdiv==0) print(i)

} # end i

R2lo <- quantile(R2.vec, probs=0.025, na.rm=T)
R2up <- quantile(R2.vec, probs=0.975, na.rm=T)
ERlo <- quantile(ER.vec, probs=0.025, na.rm=T)
ERup <- quantile(ER.vec, probs=0.975, na.rm=T)

print(c(R2lo, R2up))
print(c(ERlo, ERup))

p.NW.lo <- quantile(p.NW.vec, probs=0.025, na.rm=T)
p.NW.up <- quantile(p.NW.vec, probs=0.975, na.rm=T)
R2.NW.lo <- quantile(R2.NW.vec, probs=0.025, na.rm=T)
R2.NW.up <- quantile(R2.NW.vec, probs=0.975, na.rm=T)
slope.NW.lo <- quantile(slope.NW.vec, probs=0.025, na.rm=T)
slope.NW.up <- quantile(slope.NW.vec, probs=0.975, na.rm=T)

print(c(R2.NW.lo, R2.NW.up))
print(c(p.NW.lo, p.NW.up))
print(c(slope.NW.lo, slope.NW.up))


############################################
## regional Ricker & Gompertz logistic fits
popreg <- read.csv("popregions.csv", header=T)
head(popreg)

# EAST AND SOUTH-EASTERN ASIA
head(popreg)
ESEA <- popreg[,c(1,8)]
ESEA.Nsc <- scale(ESEA[,2], scale=T, center=F)
ESEA$Nsc <- as.numeric(ESEA.Nsc)
ESEA$r <- c(log(ESEA$Nsc[2:dim(popreg)[1]] / ESEA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(ESEA)
ESEA.use <- ESEA[which(ESEA$year >= 1962),]
plot(ESEA.use$Nsc, ESEA.use$r, xlab="scaled N", ylab="r", pch=19, ylim = c(-0.001,0.030), xlim = c(0.6, 1.5))
abline(h=0, lty=3)

# export
write.csv(ESEA.use, "ESA_use.csv", row.names=F)

# Ricker
ESEA.rick <- lm(r ~ Nsc, data=ESEA.use)
abline(ESEA.rick, lty=2, col="red")
summary(ESEA.rick)

ESEA.acf.pacf <- plot_combined_acf_pacf(na.omit(ESEA$r), max_lag = 40)
ESEA.acf.pacf

## Newey-West standard errors
ESEA.coeftest <- coeftest(ESEA.rick, vcov = NeweyWest(ESEA.rick, prewhite = FALSE), save=T)
ESEA.coeftest
ESEA.coeftest.confint <- confint(ESEA.coeftest)
ESEA.coeftest.confint

ESEA.coeftest.fitted <- attr(ESEA.coeftest, "object")[[5]]
ESEA.R2.NW <- cor(na.omit(ESEA.use$r), ESEA.coeftest.fitted)^2
ESEA.R2.NW

# Gompertz
ESEA.lNsc <- log(ESEA.use$Nsc)
ESEA.gomp <- lm(ESEA.use$r ~ ESEA.lNsc)
summary(ESEA.gomp)

ESEA.RLGL2.AIC.vec <- c(AICc(ESEA.rick), AICc(ESEA.gomp))
ESEA.RLGL2.dAIC.vec <- delta.IC(ESEA.RLGL2.AIC.vec)
ESEA.RLGL2.wAIC.vec <- weight.IC(ESEA.RLGL2.dAIC.vec)
print(ESEA.RLGL2.wAIC.vec)
ESEA.ER.RLGL2 <- ESEA.RLGL2.wAIC.vec[1]/ESEA.RLGL2.wAIC.vec[2]
ESEA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(ESEA.use$Nsc, ESEA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
ESEA.x.pred.R <- data.frame(Nsc = seq(min(ESEA.use$Nsc,na.rm=T), max(ESEA.use$Nsc,na.rm=T), length.out=100))
ESEA.pred.interval.R <- predict(ESEA.rick, newdata=ESEA.x.pred.R, interval="confidence", level=0.95)
lines(ESEA.x.pred.R$Nsc, ESEA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(ESEA.x.pred.R$Nsc, ESEA.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
ESEA.x.pred.G <- data.frame(ESEA.lNsc = seq(min(ESEA.lNsc,na.rm=T), max(ESEA.lNsc,na.rm=T), length.out=100))
ESEA.pred.interval.G <- predict(ESEA.gomp, newdata=ESEA.x.pred.G, interval="confidence", level=0.95)
lines(exp(ESEA.x.pred.G$ESEA.lNsc), ESEA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(ESEA.x.pred.G$ESEA.lNsc), ESEA.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
ESEA.R.N.descaled <- attr(ESEA.Nsc, 'scaled:scale') * ESEA.x.pred.R$Nsc
ESEA.G.N.descaled <- attr(ESEA.Nsc, 'scaled:scale') * exp(ESEA.x.pred.G$ESEA.lNsc)

## export results
ESEA.results.out <- data.frame(x.R=ESEA.R.N.descaled,
                                 R.fit=ESEA.pred.interval.R[,"fit"],
                                 R.upr=ESEA.pred.interval.R[,"upr"],
                                 R.lwr=ESEA.pred.interval.R[,"lwr"],
                                 x.G=ESEA.G.N.descaled,
                                 G.fit=ESEA.pred.interval.G[,"fit"],
                                 G.upr=ESEA.pred.interval.G[,"upr"],
                                 G.lwr=ESEA.pred.interval.G[,"lwr"])
head(ESEA.results.out)

## confidence interval for K
ESEA.negphase.st.yr <- ESEA.use$year[1]
ESEA.use.raw <- popreg[which(popreg$year >= ESEA.negphase.st.yr),]
ESEA.n.raw <- as.numeric(ESEA.use.raw$ESEA/1e3)
ESEA.r.raw <- as.numeric(ESEA.use.raw$ESEAr)
ESEA.rick.raw <- lm(ESEA.r.raw ~ ESEA.n.raw)
ESEA.confint.raw <- confint(ESEA.rick.raw)
ESEA.k.lo <- -ESEA.confint.raw[1,1]/ESEA.confint.raw[2,1]
ESEA.k.up <- -ESEA.confint.raw[1,2]/ESEA.confint.raw[2,2]
print(round(c(ESEA.k.lo, ESEA.k.up), 2))

# save to .csv
write.csv(ESEA.results.out, "ESEA_Ricker_Gompertz_fits.csv", row.names=F)


# SUB-SAHARAN AFRICA
head(popreg)
SSA <- popreg[,c(1,2)]
SSA.Nsc <- scale(SSA[,2], scale=T, center=F)
SSA$Nsc <- as.numeric(SSA.Nsc)
SSA$r <- c(log(SSA$Nsc[2:dim(popreg)[1]] / SSA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(SSA)
plot(SSA$Nsc, SSA$r, xlab="scaled N", ylab="r", pch=19)
SSA.use <- SSA[which(SSA$year >= 2008),]
plot(SSA.use$Nsc, SSA.use$r, xlab="scaled N", ylab="r", pch=19)

# export
write.csv(SSA.use, "SSA_use.csv", row.names=F)

# Ricker
SSA.rick <- lm(r ~ Nsc, data=SSA.use)
abline(SSA.rick, lty=2, col="red")
summary(SSA.rick)

SSA.acf.pacf <- plot_combined_acf_pacf(na.omit(SSA$r), max_lag = 40)
SSA.acf.pacf

## Newey-West standard errors
SSA.coeftest <- coeftest(SSA.rick, vcov = NeweyWest(SSA.rick, prewhite = FALSE), save=T)
SSA.coeftest
SSA.coeftest.confint <- confint(SSA.coeftest)
SSA.coeftest.confint
SSA.k.lo <- -SSA.coeftest.confint[1,1]/SSA.coeftest.confint[2,1]
SSA.k.up <- -SSA.coeftest.confint[1,2]/SSA.coeftest.confint[2,2]
print(round(c(SSA.k.lo, SSA.k.up), 2))

SSA.coeftest.fitted <- attr(SSA.coeftest, "object")[[5]]
SSA.R2.NW <- cor(na.omit(SSA.use$r), SSA.coeftest.fitted)^2
SSA.R2.NW

# Gompertz
SSA.lNsc <- log(SSA.use$Nsc)
SSA.gomp <- lm(SSA.use$r ~ SSA.lNsc)
summary(SSA.gomp)
SSA.RLGL2.AIC.vec <- c(AICc(SSA.rick), AICc(SSA.gomp))
SSA.RLGL2.dAIC.vec <- delta.IC(SSA.RLGL2.AIC.vec)
SSA.RLGL2.wAIC.vec <- weight.IC(SSA.RLGL2.dAIC.vec)
print(SSA.RLGL2.wAIC.vec)
SSA.ER.RLGL2 <- SSA.RLGL2.wAIC.vec[1]/SSA.RLGL2.wAIC.vec[2]
SSA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(SSA.use$Nsc, SSA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
SSA.x.pred.R <- data.frame(Nsc = seq(min(SSA.use$Nsc,na.rm=T), max(SSA.use$Nsc,na.rm=T), length.out=100))
SSA.pred.interval.R <- predict(SSA.rick, newdata=SSA.x.pred.R, interval="confidence", level=0.95)
lines(SSA.x.pred.R$Nsc, SSA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(SSA.x.pred.R$Nsc, SSA.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
SSA.x.pred.G <- data.frame(SSA.lNsc = seq(min(SSA.lNsc,na.rm=T), max(SSA.lNsc,na.rm=T), length.out=100))
SSA.pred.interval.G <- predict(SSA.gomp, newdata=SSA.x.pred.G, interval="confidence", level=0.95)
lines(exp(SSA.x.pred.G$SSA.lNsc), SSA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(SSA.x.pred.G$SSA.lNsc), SSA.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
SSA.R.N.descaled <- attr(SSA.Nsc, 'scaled:scale') * SSA.x.pred.R$Nsc
SSA.G.N.descaled <- attr(SSA.Nsc, 'scaled:scale') * exp(SSA.x.pred.G$SSA.lNsc)

## export results
SSA.results.out <- data.frame(x.R=SSA.R.N.descaled,
                              R.fit=SSA.pred.interval.R[,"fit"],
                              R.upr=SSA.pred.interval.R[,"upr"],
                              R.lwr=SSA.pred.interval.R[,"lwr"],
                              x.G=SSA.G.N.descaled,
                              G.fit=SSA.pred.interval.G[,"fit"],
                              G.upr=SSA.pred.interval.G[,"upr"],
                              G.lwr=SSA.pred.interval.G[,"lwr"])
head(SSA.results.out)

## confidence interval for K
SSA.negphase.st.yr <- SSA.use$year[1]
SSA.use.raw <- popreg[which(popreg$year >= SSA.negphase.st.yr),]
SSA.n.raw <- as.numeric(SSA.use.raw$SSA/1e3)
SSA.r.raw <- as.numeric(SSA.use.raw$SSAr)
SSA.rick.raw <- lm(SSA.r.raw ~ SSA.n.raw)
SSA.confint.raw <- confint(SSA.rick.raw)
SSA.k.lo <- -SSA.confint.raw[1,1]/SSA.confint.raw[2,1]
SSA.k.up <- -SSA.confint.raw[1,2]/SSA.confint.raw[2,2]
print(round(c(SSA.k.lo, SSA.k.up), 2))

# save to .csv
write.csv(SSA.results.out, "SSA_Ricker_Gompertz_fits.csv", row.names=F)
          
                               
# NORTH AFRICA AND WESTERN ASIA
head(popreg)
NAWA <- popreg[,c(1,4)]
NAWA.Nsc <- scale(NAWA[,2], scale=T, center=F)
NAWA$Nsc <- as.numeric(NAWA.Nsc)
NAWA$r <- c(log(NAWA$Nsc[2:dim(popreg)[1]] / NAWA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(NAWA)
plot(NAWA$Nsc, NAWA$r, xlab="scaled N", ylab="r", pch=19)
NAWA.use <- NAWA[which(NAWA$year >= 1975),]
plot(NAWA.use$Nsc, NAWA.use$r, xlab="scaled N", ylab="r", pch=19)

# export
write.csv(NAWA.use, "NAWA_use.csv", row.names=F)

# Ricker
NAWA.rick <- lm(r ~ Nsc, data=NAWA.use)
abline(NAWA.rick, lty=2, col="red")
summary(NAWA.rick)

NAWA.acf.pacf <- plot_combined_acf_pacf(na.omit(NAWA$r), max_lag = 40)
NAWA.acf.pacf

## Newey-West standard errors
NAWA.coeftest <- coeftest(NAWA.rick, vcov = NeweyWest(NAWA.rick, prewhite = FALSE), save=T)
NAWA.coeftest
NAWA.coeftest.confint <- confint(NAWA.coeftest)
NAWA.coeftest.confint

NAWA.coeftest.fitted <- attr(NAWA.coeftest, "object")[[5]]
NAWA.R2.NW <- cor(na.omit(NAWA.use$r), NAWA.coeftest.fitted)^2
NAWA.R2.NW

# Gompertz
NAWA.lNsc <- log(NAWA.use$Nsc)
NAWA.gomp <- lm(NAWA.use$r ~ NAWA.lNsc)
summary(NAWA.gomp)

NAWA.coeftest.G <- coeftest(NAWA.gomp, vcov = NeweyWest(NAWA.gomp, prewhite = FALSE), save=T)
NAWA.coeftest.G
NAWA.coeftest.confint.G <- confint(NAWA.coeftest.G)
NAWA.coeftest.confint.G

NAWA.coeftest.fitted.G <- attr(NAWA.coeftest.G, "object")[[5]]
NAWA.R2.NW.G <- cor(na.omit(NAWA.use$r), NAWA.coeftest.fitted.G)^2
NAWA.R2.NW.G

NAWA.RLGL2.AIC.vec <- c(AICc(NAWA.rick), AICc(NAWA.gomp))
NAWA.RLGL2.dAIC.vec <- delta.IC(NAWA.RLGL2.AIC.vec)
NAWA.RLGL2.wAIC.vec <- weight.IC(NAWA.RLGL2.dAIC.vec)
print(NAWA.RLGL2.wAIC.vec)
NAWA.ER.RLGL2 <- NAWA.RLGL2.wAIC.vec[1]/NAWA.RLGL2.wAIC.vec[2]
NAWA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(NAWA.use$Nsc, NAWA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
NAWA.x.pred.R <- data.frame(Nsc = seq(min(NAWA.use$Nsc,na.rm=T), max(NAWA.use$Nsc,na.rm=T), length.out=100))
NAWA.pred.interval.R <- predict(NAWA.rick, newdata=NAWA.x.pred.R, interval="confidence", level=0.95)
lines(NAWA.x.pred.R$Nsc, NAWA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(NAWA.x.pred.R$Nsc, NAWA.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
NAWA.x.pred.G <- data.frame(NAWA.lNsc = seq(min(NAWA.lNsc,na.rm=T), max(NAWA.lNsc,na.rm=T), length.out=100))
NAWA.pred.interval.G <- predict(NAWA.gomp, newdata=NAWA.x.pred.G, interval="confidence", level=0.95)
lines(exp(NAWA.x.pred.G$NAWA.lNsc), NAWA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(NAWA.x.pred.G$NAWA.lNsc), NAWA.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
NAWA.R.N.descaled <- attr(NAWA.Nsc, 'scaled:scale') * NAWA.x.pred.R$Nsc
NAWA.G.N.descaled <- attr(NAWA.Nsc, 'scaled:scale') * exp(NAWA.x.pred.G$NAWA.lNsc)

## export results
NAWA.results.out <- data.frame(x.R=NAWA.R.N.descaled,
                              R.fit=NAWA.pred.interval.R[,"fit"],
                              R.upr=NAWA.pred.interval.R[,"upr"],
                              R.lwr=NAWA.pred.interval.R[,"lwr"],
                              x.G=NAWA.G.N.descaled,
                              G.fit=NAWA.pred.interval.G[,"fit"],
                              G.upr=NAWA.pred.interval.G[,"upr"],
                              G.lwr=NAWA.pred.interval.G[,"lwr"])
head(NAWA.results.out)

## confidence interval for K
NAWA.negphase.st.yr <- NAWA.use$year[1]
NAWA.use.raw <- popreg[which(popreg$year >= NAWA.negphase.st.yr),]
NAWA.n.raw <- as.numeric(NAWA.use.raw$NAWA/1e3)
NAWA.r.raw <- as.numeric(NAWA.use.raw$NAWAr)
NAWA.rick.raw <- lm(NAWA.r.raw ~ NAWA.n.raw)
NAWA.confint.raw <- confint(NAWA.rick.raw)
NAWA.k.lo <- -NAWA.confint.raw[1,1]/NAWA.confint.raw[2,1]
NAWA.k.up <- -NAWA.confint.raw[1,2]/NAWA.confint.raw[2,2]
print(round(c(NAWA.k.lo, NAWA.k.up), 2))

# save to .csv
write.csv(NAWA.results.out, "NAWA_Ricker_Gompertz_fits.csv", row.names=F)


# CENTRAL AND SOUTH ASIA
head(popreg)
CSA <- popreg[,c(1,6)]
CSA.Nsc <- scale(CSA[,2], scale=T, center=F)
CSA$Nsc <- as.numeric(CSA.Nsc)
CSA$r <- c(log(CSA$Nsc[2:dim(popreg)[1]] / CSA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(CSA)
plot(CSA$Nsc, CSA$r, xlab="scaled N", ylab="r", pch=19)
CSA.use <- CSA[which(CSA$year >= 1983),]
plot(CSA.use$Nsc, CSA.use$r, xlab="scaled N", ylab="r", pch=19)

# export
write.csv(CSA.use, "CSA_use.csv", row.names=F)

# Ricker
CSA.rick <- lm(r ~ Nsc, data=CSA.use)
abline(CSA.rick, lty=2, col="red")
summary(CSA.rick)

CSA.acf.pacf <- plot_combined_acf_pacf(na.omit(CSA$r), max_lag = 40)
CSA.acf.pacf

## Newey-West standard errors
CSA.coeftest <- coeftest(CSA.rick, vcov = NeweyWest(CSA.rick, prewhite = FALSE), save=T)
CSA.coeftest
CSA.coeftest.confint <- confint(CSA.coeftest)
CSA.coeftest.confint

CSA.coeftest.fitted <- attr(CSA.coeftest, "object")[[5]]
CSA.R2.NW <- cor(na.omit(CSA.use$r), CSA.coeftest.fitted)^2
CSA.R2.NW

# Gompertz
CSA.lNsc <- log(CSA.use$Nsc)
CSA.gomp <- lm(CSA.use$r ~ CSA.lNsc)
summary(CSA.gomp)

CSA.RLGL2.AIC.vec <- c(AICc(CSA.rick), AICc(CSA.gomp))
CSA.RLGL2.dAIC.vec <- delta.IC(CSA.RLGL2.AIC.vec)
CSA.RLGL2.wAIC.vec <- weight.IC(CSA.RLGL2.dAIC.vec)
print(CSA.RLGL2.wAIC.vec)
CSA.ER.RLGL2 <- CSA.RLGL2.wAIC.vec[1]/CSA.RLGL2.wAIC.vec[2]
CSA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(CSA.use$Nsc, CSA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
CSA.x.pred.R <- data.frame(Nsc = seq(min(CSA.use$Nsc,na.rm=T), max(CSA.use$Nsc,na.rm=T), length.out=100))
CSA.pred.interval.R <- predict(CSA.rick, newdata=CSA.x.pred.R, interval="confidence", level=0.95)
lines(CSA.x.pred.R$Nsc, CSA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(CSA.x.pred.R$Nsc, CSA.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
CSA.x.pred.G <- data.frame(CSA.lNsc = seq(min(CSA.lNsc,na.rm=T), max(CSA.lNsc,na.rm=T), length.out=100))
CSA.pred.interval.G <- predict(CSA.gomp, newdata=CSA.x.pred.G, interval="confidence", level=0.95)
lines(exp(CSA.x.pred.G$CSA.lNsc), CSA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(CSA.x.pred.G$CSA.lNsc), CSA.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
CSA.R.N.descaled <- attr(CSA.Nsc, 'scaled:scale') * CSA.x.pred.R$Nsc
CSA.G.N.descaled <- attr(CSA.Nsc, 'scaled:scale') * exp(CSA.x.pred.G$CSA.lNsc)

## export results
CSA.results.out <- data.frame(x.R=CSA.R.N.descaled,
                              R.fit=CSA.pred.interval.R[,"fit"],
                              R.upr=CSA.pred.interval.R[,"upr"],
                              R.lwr=CSA.pred.interval.R[,"lwr"],
                              x.G=CSA.G.N.descaled,
                              G.fit=CSA.pred.interval.G[,"fit"],
                              G.upr=CSA.pred.interval.G[,"upr"],
                              G.lwr=CSA.pred.interval.G[,"lwr"])
head(CSA.results.out)

## confidence interval for K
CSA.negphase.st.yr <- CSA.use$year[1]
CSA.use.raw <- popreg[which(popreg$year >= CSA.negphase.st.yr),]
CSA.n.raw <- as.numeric(CSA.use.raw$CSA/1e3)
CSA.r.raw <- as.numeric(CSA.use.raw$CSAr)
CSA.rick.raw <- lm(CSA.r.raw ~ CSA.n.raw)
CSA.confint.raw <- confint(CSA.rick.raw)
CSA.k.lo <- -CSA.confint.raw[1,1]/CSA.confint.raw[2,1]
CSA.k.up <- -CSA.confint.raw[1,2]/CSA.confint.raw[2,2]
print(round(c(CSA.k.lo, CSA.k.up), 2))

# save to .csv
write.csv(CSA.results.out, "CSA_Ricker_Gompertz_fits.csv", row.names=F)



# LATIN AMERICA AND CARIBBEAN
head(popreg)
LAC <- popreg[,c(1,10)]
LAC.Nsc <- scale(LAC[,2], scale=T, center=F)
LAC$Nsc <- as.numeric(LAC.Nsc)
LAC$r <- c(log(LAC$Nsc[2:dim(popreg)[1]] / LAC$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(LAC)
plot(LAC$Nsc, LAC$r, xlab="scaled N", ylab="r", pch=19)
LAC.use <- LAC[which(LAC$year >= 1961),]
plot(LAC.use$Nsc, LAC.use$r, xlab="scaled N", ylab="r", pch=19)

# export
write.csv(LAC.use, "LAC_use.csv", row.names=F)

# Ricker
LAC.rick <- lm(r ~ Nsc, data=LAC.use)
abline(LAC.rick, lty=2, col="red")
summary(LAC.rick)

LAC.acf.pacf <- plot_combined_acf_pacf(na.omit(LAC$r), max_lag = 40)
LAC.acf.pacf

## Newey-West standard errors
LAC.coeftest <- coeftest(LAC.rick, vcov = NeweyWest(LAC.rick, prewhite = FALSE), save=T)
LAC.coeftest
LAC.coeftest.confint <- confint(LAC.coeftest)
LAC.coeftest.confint

LAC.coeftest.fitted <- attr(LAC.coeftest, "object")[[5]]
LAC.R2.NW <- cor(na.omit(LAC.use$r), LAC.coeftest.fitted)^2
LAC.R2.NW

# Gompertz
LAC.lNsc <- log(LAC.use$Nsc)
LAC.gomp <- lm(LAC.use$r ~ LAC.lNsc)
summary(LAC.gomp)

LAC.RLGL2.AIC.vec <- c(AICc(LAC.rick), AICc(LAC.gomp))
LAC.RLGL2.dAIC.vec <- delta.IC(LAC.RLGL2.AIC.vec)
LAC.RLGL2.wAIC.vec <- weight.IC(LAC.RLGL2.dAIC.vec)
print(LAC.RLGL2.wAIC.vec)
LAC.ER.RLGL2 <- LAC.RLGL2.wAIC.vec[1]/LAC.RLGL2.wAIC.vec[2]
LAC.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(LAC.use$Nsc, LAC.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
LAC.x.pred.R <- data.frame(Nsc = seq(min(LAC.use$Nsc,na.rm=T), max(LAC.use$Nsc,na.rm=T), length.out=100))
LAC.pred.interval.R <- predict(LAC.rick, newdata=LAC.x.pred.R, interval="confidence", level=0.95)
lines(LAC.x.pred.R$Nsc, LAC.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(LAC.x.pred.R$Nsc, LAC.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
LAC.x.pred.G <- data.frame(LAC.lNsc = seq(min(LAC.lNsc,na.rm=T), max(LAC.lNsc,na.rm=T), length.out=100))
LAC.pred.interval.G <- predict(LAC.gomp, newdata=LAC.x.pred.G, interval="confidence", level=0.95)
lines(exp(LAC.x.pred.G$LAC.lNsc), LAC.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(LAC.x.pred.G$LAC.lNsc), LAC.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
LAC.R.N.descaled <- attr(LAC.Nsc, 'scaled:scale') * LAC.x.pred.R$Nsc
LAC.G.N.descaled <- attr(LAC.Nsc, 'scaled:scale') * exp(LAC.x.pred.G$LAC.lNsc)

## export results
LAC.results.out <- data.frame(x.R=LAC.R.N.descaled,
                              R.fit=LAC.pred.interval.R[,"fit"],
                              R.upr=LAC.pred.interval.R[,"upr"],
                              R.lwr=LAC.pred.interval.R[,"lwr"],
                              x.G=LAC.G.N.descaled,
                              G.fit=LAC.pred.interval.G[,"fit"],
                              G.upr=LAC.pred.interval.G[,"upr"],
                              G.lwr=LAC.pred.interval.G[,"lwr"])
head(LAC.results.out)

## confidence interval for K
LAC.negphase.st.yr <- LAC.use$year[1]
LAC.use.raw <- popreg[which(popreg$year >= LAC.negphase.st.yr),]
LAC.n.raw <- as.numeric(LAC.use.raw$LAC/1e3)
LAC.r.raw <- as.numeric(LAC.use.raw$LACr)
LAC.rick.raw <- lm(LAC.r.raw ~ LAC.n.raw)
LAC.confint.raw <- confint(LAC.rick.raw)
LAC.k.lo <- -LAC.confint.raw[1,1]/LAC.confint.raw[2,1]
LAC.k.up <- -LAC.confint.raw[1,2]/LAC.confint.raw[2,2]
print(round(c(LAC.k.lo, LAC.k.up), 2))

# save to .csv
write.csv(LAC.results.out, "LAC_Ricker_Gompertz_fits.csv", row.names=F)



# OCEANIA
head(popreg)
OC <- popreg[,c(1,12)]
OC.Nsc <- scale(OC[,2], scale=T, center=F)
OC$Nsc <- as.numeric(OC.Nsc)
OC$r <- c(log(OC$Nsc[2:dim(popreg)[1]] / OC$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(OC)
plot(OC$Nsc, OC$r, xlab="scaled N", ylab="r", pch=19)
OC.use <- OC[which(OC$year >= 2016),]
plot(OC.use$Nsc, OC.use$r, xlab="scaled N", ylab="r", pch=19)

# export
write.csv(OC.use, "OC_use.csv", row.names=F)

# Ricker
OC.rick <- lm(r ~ Nsc, data=OC.use)
abline(OC.rick, lty=2, col="red")
summary(OC.rick)
OC.acf.pacf <- plot_combined_acf_pacf(na.omit(OC$r), max_lag = 40)
OC.acf.pacf

## Newey-West standard errors
OC.coeftest <- coeftest(OC.rick, vcov = NeweyWest(OC.rick, prewhite = FALSE), save=T)
OC.coeftest
OC.coeftest.confint <- confint(OC.coeftest)
OC.coeftest.confint

OC.coeftest.fitted <- attr(OC.coeftest, "object")[[5]]
OC.R2.NW <- cor(na.omit(OC.use$r), OC.coeftest.fitted)^2
OC.R2.NW

# Gompertz
OC.lNsc <- log(OC.use$Nsc)
OC.gomp <- lm(OC.use$r ~ OC.lNsc)
summary(OC.gomp)

OC.RLGL2.AIC.vec <- c(AICc(OC.rick), AICc(OC.gomp))
OC.RLGL2.dAIC.vec <- delta.IC(OC.RLGL2.AIC.vec)
OC.RLGL2.wAIC.vec <- weight.IC(OC.RLGL2.dAIC.vec)
print(OC.RLGL2.wAIC.vec)
OC.ER.RLGL2 <- OC.RLGL2.wAIC.vec[1]/OC.RLGL2.wAIC.vec[2]
OC.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(OC.use$Nsc, OC.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
OC.x.pred.R <- data.frame(Nsc = seq(min(OC.use$Nsc,na.rm=T), max(OC.use$Nsc,na.rm=T), length.out=100))
OC.pred.interval.R <- predict(OC.rick, newdata=OC.x.pred.R, interval="confidence", level=0.95)
lines(OC.x.pred.R$Nsc, OC.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(OC.x.pred.R$Nsc, OC.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
OC.x.pred.G <- data.frame(OC.lNsc = seq(min(OC.lNsc,na.rm=T), max(OC.lNsc,na.rm=T), length.out=100))
OC.pred.interval.G <- predict(OC.gomp, newdata=OC.x.pred.G, interval="confidence", level=0.95)
lines(exp(OC.x.pred.G$OC.lNsc), OC.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(OC.x.pred.G$OC.lNsc), OC.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
OC.R.N.descaled <- attr(OC.Nsc, 'scaled:scale') * OC.x.pred.R$Nsc
OC.G.N.descaled <- attr(OC.Nsc, 'scaled:scale') * exp(OC.x.pred.G$OC.lNsc)

## export results
OC.results.out <- data.frame(x.R=OC.R.N.descaled,
                              R.fit=OC.pred.interval.R[,"fit"],
                              R.upr=OC.pred.interval.R[,"upr"],
                              R.lwr=OC.pred.interval.R[,"lwr"],
                              x.G=OC.G.N.descaled,
                              G.fit=OC.pred.interval.G[,"fit"],
                              G.upr=OC.pred.interval.G[,"upr"],
                              G.lwr=OC.pred.interval.G[,"lwr"])
head(OC.results.out)

## confidence interval for K
OC.negphase.st.yr <- OC.use$year[1]
OC.use.raw <- popreg[which(popreg$year >= OC.negphase.st.yr),]
OC.n.raw <- as.numeric(OC.use.raw$OC/1e3)
OC.r.raw <- as.numeric(OC.use.raw$OCr)
OC.rick.raw <- lm(OC.r.raw ~ OC.n.raw)
OC.confint.raw <- confint(OC.rick.raw)
OC.k.lo <- -OC.confint.raw[1,1]/OC.confint.raw[2,1]
OC.k.up <- -OC.confint.raw[1,2]/OC.confint.raw[2,2]
print(round(c(OC.k.lo, OC.k.up), 2))


# save to .csv
write.csv(OC.results.out, "OC_Ricker_Gompertz_fits.csv", row.names=F)


# OCEANIA excluding Australian and NZ
head(popreg)
OC.noANZ <- popreg[,c(1,16)]
OC.noANZ.Nsc <- scale(OC.noANZ[,2], scale=T, center=F)
OC.noANZ$Nsc <- as.numeric(OC.noANZ.Nsc)
OC.noANZ$r <- c(log(OC.noANZ$Nsc[2:dim(popreg)[1]] / OC.noANZ$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(OC.noANZ)
plot(OC.noANZ$Nsc, OC.noANZ$r, xlab="scaled N", ylab="r", pch=19)
OC.noANZ.use <- OC.noANZ[which(OC.noANZ$year >= 1996),]
plot(OC.noANZ.use$Nsc, OC.noANZ.use$r, xlab="scaled N", ylab="r", pch=19)

# export
write.csv(OC.noANZ.use, "OCnoANZ_use.csv", row.names=F)

# Ricker
OC.noANZ.rick <- lm(r ~ Nsc, data=OC.noANZ.use)
abline(OC.noANZ.rick, lty=2, col="red")
summary(OC.noANZ.rick)
OC.noANZ.acf.pacf <- plot_combined_acf_pacf(na.omit(OC.noANZ$r), max_lag = 40)
OC.noANZ.acf.pacf

## Newey-West standard errors
OC.noANZ.coeftest <- coeftest(OC.noANZ.rick, vcov = NeweyWest(OC.noANZ.rick, prewhite = FALSE), save=T)
OC.noANZ.coeftest
OC.noANZ.coeftest.confint <- confint(OC.noANZ.coeftest)
OC.noANZ.coeftest.confint

OC.noANZ.coeftest.fitted <- attr(OC.noANZ.coeftest, "object")[[5]]
OC.noANZ.R2.NW <- cor(na.omit(OC.noANZ.use$r), OC.noANZ.coeftest.fitted)^2
OC.noANZ.R2.NW

# Gompertz
OC.noANZ.lNsc <- log(OC.noANZ.use$Nsc)
OC.noANZ.gomp <- lm(OC.noANZ.use$r ~ OC.noANZ.lNsc)
summary(OC.noANZ.gomp)

OC.noANZ.RLGL2.AIC.vec <- c(AICc(OC.noANZ.rick), AICc(OC.noANZ.gomp))
OC.noANZ.RLGL2.dAIC.vec <- delta.IC(OC.noANZ.RLGL2.AIC.vec)
OC.noANZ.RLGL2.wAIC.vec <- weight.IC(OC.noANZ.RLGL2.dAIC.vec)
print(OC.noANZ.RLGL2.wAIC.vec)
OC.noANZ.ER.RLGL2 <- OC.noANZ.RLGL2.wAIC.vec[1]/OC.noANZ.RLGL2.wAIC.vec[2]
OC.noANZ.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(OC.noANZ.use$Nsc, OC.noANZ.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
OC.noANZ.x.pred.R <- data.frame(Nsc = seq(min(OC.noANZ.use$Nsc,na.rm=T), max(OC.noANZ.use$Nsc,na.rm=T), length.out=100))
OC.noANZ.pred.interval.R <- predict(OC.noANZ.rick, newdata=OC.noANZ.x.pred.R, interval="confidence", level=0.95)
lines(OC.noANZ.x.pred.R$Nsc, OC.noANZ.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(OC.noANZ.x.pred.R$Nsc, OC.noANZ.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
OC.noANZ.x.pred.G <- data.frame(OC.noANZ.lNsc = seq(min(OC.noANZ.lNsc,na.rm=T), max(OC.noANZ.lNsc,na.rm=T), length.out=100))
OC.noANZ.pred.interval.G <- predict(OC.noANZ.gomp, newdata=OC.noANZ.x.pred.G, interval="confidence", level=0.95)
lines(exp(OC.noANZ.x.pred.G$OC.noANZ.lNsc), OC.noANZ.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(OC.noANZ.x.pred.G$OC.noANZ.lNsc), OC.noANZ.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
OC.noANZ.R.N.descaled <- attr(OC.noANZ.Nsc, 'scaled:scale') * OC.noANZ.x.pred.R$Nsc
OC.noANZ.G.N.descaled <- attr(OC.noANZ.Nsc, 'scaled:scale') * exp(OC.noANZ.x.pred.G$OC.noANZ.lNsc)

## export results
OC.noANZ.results.out <- data.frame(x.R=OC.noANZ.R.N.descaled,
                              R.fit=OC.noANZ.pred.interval.R[,"fit"],
                              R.upr=OC.noANZ.pred.interval.R[,"upr"],
                              R.lwr=OC.noANZ.pred.interval.R[,"lwr"],
                              x.G=OC.noANZ.G.N.descaled,
                              G.fit=OC.noANZ.pred.interval.G[,"fit"],
                              G.upr=OC.noANZ.pred.interval.G[,"upr"],
                              G.lwr=OC.noANZ.pred.interval.G[,"lwr"])
head(OC.noANZ.results.out)

## confidence interval for K
OC.noANZ.negphase.st.yr <- OC.noANZ.use$year[1]
OC.noANZ.use.raw <- popreg[which(popreg$year >= OC.noANZ.negphase.st.yr),]
OC.noANZ.n.raw <- as.numeric(OC.noANZ.use.raw$OCnANZ/1e3)
OC.noANZ.r.raw <- as.numeric(OC.noANZ.use.raw$OCnANZr)
OC.noANZ.rick.raw <- lm(OC.noANZ.r.raw ~ OC.noANZ.n.raw)
OC.noANZ.confint.raw <- confint(OC.noANZ.rick.raw)
OC.noANZ.k.lo <- -OC.noANZ.confint.raw[1,1]/OC.noANZ.confint.raw[2,1]
OC.noANZ.k.up <- -OC.noANZ.confint.raw[1,2]/OC.noANZ.confint.raw[2,2]
print(round(c(OC.noANZ.k.lo, OC.noANZ.k.up), 2))

# save to .csv
write.csv(OC.noANZ.results.out, "OC_noANZ_Ricker_Gompertz_fits.csv", row.names=F)




# EUROPE AND NORTH AMERICA
head(popreg)
EUNA <- popreg[,c(1,14)]
EUNA.Nsc <- scale(EUNA[,2], scale=T, center=F)
EUNA$Nsc <- as.numeric(EUNA.Nsc)
EUNA$r <- c(log(EUNA$Nsc[2:dim(popreg)[1]] / EUNA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(EUNA)
plot(EUNA$Nsc, EUNA$r, xlab="scaled N", ylab="r", pch=19)
EUNA.use <- EUNA[which(EUNA$year >= 1957),]
plot(EUNA.use$Nsc, EUNA.use$r, xlab="scaled N", ylab="r", pch=19)

# export
write.csv(EUNA.use, "EUNA_use.csv", row.names=F)

# Ricker
EUNA.rick <- lm(r ~ Nsc, data=EUNA.use)
abline(EUNA.rick, lty=2, col="red")
summary(EUNA.rick)

EUNA.acf.pacf <- plot_combined_acf_pacf(na.omit(EUNA$r), max_lag = 40)
EUNA.acf.pacf

## Newey-West standard errors
EUNA.coeftest <- coeftest(EUNA.rick, vcov = NeweyWest(EUNA.rick, prewhite = FALSE), save=T)
EUNA.coeftest
EUNA.coeftest.confint <- confint(EUNA.coeftest)
EUNA.coeftest.confint

EUNA.coeftest.fitted <- attr(EUNA.coeftest, "object")[[5]]
EUNA.R2.NW <- cor(na.omit(EUNA.use$r), EUNA.coeftest.fitted)^2
EUNA.R2.NW

# Gompertz
EUNA.lNsc <- log(EUNA.use$Nsc)
EUNA.gomp <- lm(EUNA.use$r ~ EUNA.lNsc)
summary(EUNA.gomp)

EUNA.RLGL2.AIC.vec <- c(AICc(EUNA.rick), AICc(EUNA.gomp))
EUNA.RLGL2.dAIC.vec <- delta.IC(EUNA.RLGL2.AIC.vec)
EUNA.RLGL2.wAIC.vec <- weight.IC(EUNA.RLGL2.dAIC.vec)
print(EUNA.RLGL2.wAIC.vec)
EUNA.ER.RLGL2 <- EUNA.RLGL2.wAIC.vec[1]/EUNA.RLGL2.wAIC.vec[2]
EUNA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(EUNA.use$Nsc, EUNA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
EUNA.x.pred.R <- data.frame(Nsc = seq(min(EUNA.use$Nsc,na.rm=T), max(EUNA.use$Nsc,na.rm=T), length.out=100))
EUNA.pred.interval.R <- predict(EUNA.rick, newdata=EUNA.x.pred.R, interval="confidence", level=0.95)
lines(EUNA.x.pred.R$Nsc, EUNA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(EUNA.x.pred.R$Nsc, EUNA.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
EUNA.x.pred.G <- data.frame(EUNA.lNsc = seq(min(EUNA.lNsc,na.rm=T), max(EUNA.lNsc,na.rm=T), length.out=100))
EUNA.pred.interval.G <- predict(EUNA.gomp, newdata=EUNA.x.pred.G, interval="confidence", level=0.95)
lines(exp(EUNA.x.pred.G$EUNA.lNsc), EUNA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(EUNA.x.pred.G$EUNA.lNsc), EUNA.pred.interval.G[,"upr"], col="green", lty=2) 

# unscale Nsc
EUNA.R.N.descaled <- attr(EUNA.Nsc, 'scaled:scale') * EUNA.x.pred.R$Nsc
EUNA.G.N.descaled <- attr(EUNA.Nsc, 'scaled:scale') * exp(EUNA.x.pred.G$EUNA.lNsc)

## export results
EUNA.results.out <- data.frame(x.R=EUNA.R.N.descaled,
                              R.fit=EUNA.pred.interval.R[,"fit"],
                              R.upr=EUNA.pred.interval.R[,"upr"],
                              R.lwr=EUNA.pred.interval.R[,"lwr"],
                              x.G=EUNA.G.N.descaled,
                              G.fit=EUNA.pred.interval.G[,"fit"],
                              G.upr=EUNA.pred.interval.G[,"upr"],
                              G.lwr=EUNA.pred.interval.G[,"lwr"])
head(EUNA.results.out)

## confidence interval for K
EUNA.negphase.st.yr <- EUNA.use$year[1]
EUNA.use.raw <- popreg[which(popreg$year >= EUNA.negphase.st.yr),]
EUNA.n.raw <- as.numeric(EUNA.use.raw$EUNA/1e3)
EUNA.r.raw <- as.numeric(EUNA.use.raw$EUNAr)
EUNA.rick.raw <- lm(EUNA.r.raw ~ EUNA.n.raw)
EUNA.confint.raw <- confint(EUNA.rick.raw)
EUNA.k.lo <- -EUNA.confint.raw[1,1]/EUNA.confint.raw[2,1]
EUNA.k.up <- -EUNA.confint.raw[1,2]/EUNA.confint.raw[2,2]
print(round(c(EUNA.k.lo, EUNA.k.up), 2))

# save to .csv
write.csv(EUNA.results.out, "EUNA_Ricker_Gompertz_fits.csv", row.names=F)


## sum regional K results
regional.K.results <- data.frame(Region=c("SSA","ESEA","NAWA","EUNA","LAC","OC","OCnANZ","CAS"),
                                   K.lo=c(SSA.k.lo, ESEA.k.lo, NAWA.k.lo, EUNA.k.lo, LAC.k.lo, OC.k.lo, CSA.k.lo),
                                   K.up=c(SSA.k.up, ESEA.k.up, NAWA.k.up, EUNA.k.up, LAC.k.up, OC.k.up, CSA.k.up))
regional.K.results

# sum columns
colSums(regional.K.results[,2:3])
reg.sum.med <- mean(colSums(regional.K.results[,2:3]))
reg.sum.med
glob.sum.med <- mean(c(k3.lo,k3.up))/10^9
glob.sum.med

# % diff
(reg.sum.med - glob.sum.med) / glob.sum.med * 100


## import regional demography time series
demog <- read.csv("demog.csv", header=T)
head(demog)

## regions vector
regs.vec <- attr(table(demog$region), 'names')
regs.vec

## year when tfr reached ≤ 2.1
decl.phas.st.yr.vec <- c(CSA.negphase.st.yr,ESEA.negphase.st.yr,EUNA.negphase.st.yr,
                         LAC.negphase.st.yr,NAWA.negphase.st.yr,OC.negphase.st.yr, 
                         OC.noANZ.negphase.st.yr,SSA.negphase.st.yr)
names(decl.phas.st.yr.vec) <- regs.vec
decl.phas.st.yr.vec

## average tfr by region
tfr.avg.vec <- c()
for (i in 1:length(regs.vec)) {
  reg.i <- regs.vec[i]
  demog.reg.i <- demog[which(demog$region == reg.i),]
  tfr.avg.i <- mean(demog.reg.i$tfr, na.rm=T)
  tfr.avg.vec <- c(tfr.avg.vec, tfr.avg.i)
}
tfr.avg.vec

## average tfr during decline phase per region
tfr.declphase.avg.vec <- c()
for (i in 1:length(regs.vec)) {
  reg.i <- regs.vec[i]
  demog.reg.i <- demog[which(demog$region == reg.i),]
  decl.phas.st.yr.i <- decl.phas.st.yr.vec[i]
  demog.reg.i.declphase <- demog.reg.i[which(demog.reg.i$year >= decl.phas.st.yr.i),]
  tfr.declphase.avg.i <- mean(demog.reg.i.declphase$tfr, na.rm=T)
  tfr.declphase.avg.vec <- c(tfr.declphase.avg.vec, tfr.declphase.avg.i)
}
tfr.declphase.avg.vec

# combine into data frame
demog.reg.summary <- data.frame(region=regs.vec, tfr.avg=tfr.avg.vec,
                                decl_st_yr=decl.phas.st.yr.vec,
                                tfr.declphase.avg=tfr.declphase.avg.vec)
demog.reg.summary

## per-capita PPP-GDP average by region
# country & regions classification
cont.cntry <- read.csv("continent.country2.csv", header=T)
head(cont.cntry)

## import pcpppgdp by country
pcpppgdp <- read.csv("gdppcPPP.csv", header=T)
head(pcpppgdp)

## merge cont.cntry and pcpppgdp
pcpppgdp.reg <- merge(pcpppgdp, cont.cntry, by="cntry.code")
head(pcpppgdp.reg)

table(pcpppgdp.reg$region)
table(demog$region)

pcpppgdp.reg$country[which(pcpppgdp.reg$region == "AUS")]

## reclassify region names to match demog.reg.summary
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "CASIA", "CSA", NA)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "EASIA", "ESEA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "SEA", "ESEA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "SASIA", "CSA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "CAM", "LAC", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "CAR", "LAC", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "SA", "LAC", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "EAFR", "SSA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "MAFR", "SSA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "WAFR", "SSA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "SAFR", "SSA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "EU", "EUNA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "NAM", "EUNA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "NEU", "EUNA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "ME", "NAWA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "NAFR", "NAWA", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "POLY", "OC", pcpppgdp.reg$region2)
pcpppgdp.reg$region2 <- ifelse(pcpppgdp.reg$region == "AUS", "OC", pcpppgdp.reg$region2)

## average pcpppgdp by country 1990-2024
pcpppgdp.reg$pcpppgdp <- rowMeans(pcpppgdp.reg[,2:36], na.rm=T)
head(pcpppgdp.reg)

## average pcpppgdp by region
pcpppgdp.reg.avg <- aggregate(pcpppgdp.reg$pcpppgdp, by=list(region=pcpppgdp.reg$region2), FUN=mean, na.rm=T)
colnames(pcpppgdp.reg.avg)[2] <- "pcpppgdp.avg"
pcpppgdp.reg.avg

## remove Australia & New Zealand and recalculate for Oceania
pcpppgdp.reg.noANZ <- pcpppgdp.reg[which(pcpppgdp.reg$cntry.code != "AUS" & pcpppgdp.reg$cntry.code != "NZ"),]
pcpppgdp.reg.noANZ.avg <- aggregate(pcpppgdp.reg.noANZ$pcpppgdp, by=list(region=pcpppgdp.reg.noANZ$region2), FUN=mean, na.rm=T)
colnames(pcpppgdp.reg.noANZ.avg)[2] <- "pcpppgdp.avg"
pcpppgdp.reg.noANZ.avg

## rate of fertility change
head(demog)
first.neg.tfr_change.yr.vec <- diff_mean_tfr_change_bef_aft <- rep(NA, length(regs.vec))
for (r in 1:length(regs.vec)) {
  reg.i <- regs.vec[r]
  demog.reg.i <- demog[which(demog$region == reg.i),]
  tfr_change <- log(demog.reg.i$tfr[2:length(demog.reg.i$tfr)] / demog.reg.i$tfr[1:(length(demog.reg.i$tfr)-1)])
  plot(demog.reg.i$year[2:length(demog.reg.i$year)], tfr_change, type="l", xlab="year", ylab="log(tfr change)",
       main=paste("region:", reg.i))
  first.neg.tfr_change.yr.vec[r] <- demog.reg.i$year[which(tfr_change < 0)[1] + 1]
  st_yr_bef <- which(demog.reg.i$year == (decl.phas.st.yr.vec[r] - 4))
  end_yr_bef <- which(demog.reg.i$year == decl.phas.st.yr.vec[r])  
  mean_tfr_change_bef <- mean(tfr_change[st_yr_bef:end_yr_bef], na.rm=T)
  end_yr_aft <- which(demog.reg.i$year == (decl.phas.st.yr.vec[r] + 4))
  mean_tfr_change_aft <- mean(tfr_change[end_yr_bef:end_yr_aft], na.rm=T)
  diff_mean_tfr_change_bef_aft[r] <- mean_tfr_change_aft - mean_tfr_change_bef
}
names(first.neg.tfr_change.yr.vec) <- regs.vec
first.neg.tfr_change.yr.vec
names(diff_mean_tfr_change_bef_aft) <- regs.vec
diff_mean_tfr_change_bef_aft

## top 10 countries with highest fertility
# Niger
NER <- read.csv("NER.csv", header=T)
head(NER)
NER.Nsc <- scale(NER[,2], scale=T, center=F)
NER$Nsc <- as.numeric(NER.Nsc)
NER$r <- c(log(NER$Nsc[2:dim(popreg)[1]] / NER$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(NER)
plot(NER$Nsc, NER$r, xlab="N", ylab="r", pch=19)
NER.use <- NER[which(NER$year >= 2012),]
plot(NER.use$Nsc, NER.use$r, xlab="N", ylab="r", pch=19)

# Ricker
NER.rick <- lm(r ~ Nsc, data=NER.use)
abline(NER.rick, lty=2, col="red")
summary(NER.rick)

NER.acf.pacf <- plot_combined_acf_pacf(na.omit(NER$r), max_lag = 40)
NER.acf.pacf

## Newey-West standard errors
NER.coeftest <- coeftest(NER.rick, vcov = NeweyWest(NER.rick, prewhite = FALSE), save=T)
NER.coeftest
NER.coeftest.confint <- confint(NER.coeftest)
NER.coeftest.confint

NER.coeftest.fitted <- attr(NER.coeftest, "object")[[5]]
NER.R2.NW <- cor(na.omit(NER.use$r), NER.coeftest.fitted)^2
NER.R2.NW

# Gompertz
NER.lNsc <- log(NER.use$Nsc)
NER.gomp <- lm(NER.use$r ~ NER.lNsc)
summary(NER.gomp)

NER.RLGL2.AIC.vec <- c(AICc(NER.rick), AICc(NER.gomp))
NER.RLGL2.dAIC.vec <- delta.IC(NER.RLGL2.AIC.vec)
NER.RLGL2.wAIC.vec <- weight.IC(NER.RLGL2.dAIC.vec)
print(NER.RLGL2.wAIC.vec)
NER.ER.RLGL2 <- NER.RLGL2.wAIC.vec[1]/NER.RLGL2.wAIC.vec[2]
NER.ER.RLGL2
1/NER.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(NER.use$Nsc, NER.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
NER.x.pred.R <- data.frame(Nsc = seq(min(NER.use$Nsc,na.rm=T), max(NER.use$Nsc,na.rm=T), length.out=100))
NER.pred.interval.R <- predict(NER.rick, newdata=NER.x.pred.R, interval="confidence", level=0.95)
lines(NER.x.pred.R$Nsc, NER.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(NER.x.pred.R$Nsc, NER.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
NER.x.pred.G <- data.frame(NER.lNsc = seq(min(NER.lNsc,na.rm=T), max(NER.lNsc,na.rm=T), length.out=100))
NER.pred.interval.G <- predict(NER.gomp, newdata=NER.x.pred.G, interval="confidence", level=0.95)
lines(exp(NER.x.pred.G$NER.lNsc), NER.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(NER.x.pred.G$NER.lNsc), NER.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
NER.R.N.descaled <- attr(NER.Nsc, 'scaled:scale') * NER.x.pred.R$Nsc
NER.G.N.descaled <- attr(NER.Nsc, 'scaled:scale') * exp(NER.x.pred.G$NER.lNsc)

## export results
NER.results.out <- data.frame(x.R=NER.R.N.descaled,
                              R.fit=NER.pred.interval.R[,"fit"],
                              R.upr=NER.pred.interval.R[,"upr"],
                              R.lwr=NER.pred.interval.R[,"lwr"],
                              x.G=NER.G.N.descaled,
                              G.fit=NER.pred.interval.G[,"fit"],
                              G.upr=NER.pred.interval.G[,"upr"],
                              G.lwr=NER.pred.interval.G[,"lwr"])
head(NER.results.out)

# save to .csv
write.csv(NER.results.out, "NER_Ricker_Gompertz_fits.csv", row.names=F)


# Democratic Republic of Congo
COD <- read.csv("COD.csv", header=T)
head(COD)
COD.Nsc <- scale(COD[,2], scale=T, center=F)
COD$Nsc <- as.numeric(COD.Nsc)
COD$r <- c(log(COD$Nsc[2:dim(popreg)[1]] / COD$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(COD)
plot(COD$Nsc, COD$r, xlab="N", ylab="r", pch=19)
COD.use <- COD[which(COD$year >= 2013),]
plot(COD.use$Nsc, COD.use$r, xlab="N", ylab="r", pch=19)                                         

# Ricker
COD.rick <- lm(r ~ Nsc, data=COD.use)
abline(COD.rick, lty=2, col="red")
summary(COD.rick)

COD.acf.pacf <- plot_combined_acf_pacf(na.omit(COD$r), max_lag = 40)
COD.acf.pacf

## Newey-West standard errors
COD.coeftest <- coeftest(COD.rick, vcov = NeweyWest(COD.rick, prewhite = FALSE), save=T)
COD.coeftest
COD.coeftest.confint <- confint(COD.coeftest)
COD.coeftest.confint

COD.coeftest.fitted <- attr(COD.coeftest, "object")[[5]]
COD.R2.NW <- cor(na.omit(COD.use$r), COD.coeftest.fitted)^2
COD.R2.NW

# Gompertz
COD.lNsc <- log(COD.use$Nsc)
COD.gomp <- lm(COD.use$r ~ COD.lNsc)
summary(COD.gomp)

COD.RLGL2.AIC.vec <- c(AICc(COD.rick), AICc(COD.gomp))
COD.RLGL2.dAIC.vec <- delta.IC(COD.RLGL2.AIC.vec)
COD.RLGL2.wAIC.vec <- weight.IC(COD.RLGL2.dAIC.vec)
print(COD.RLGL2.wAIC.vec)
COD.ER.RLGL2 <- COD.RLGL2.wAIC.vec[1]/COD.RLGL2.wAIC.vec[2]
COD.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(COD.use$Nsc, COD.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
COD.x.pred.R <- data.frame(Nsc = seq(min(COD.use$Nsc,na.rm=T), max(COD.use$Nsc,na.rm=T), length.out=100))
COD.pred.interval.R <- predict(COD.rick, newdata=COD.x.pred.R, interval="confidence", level=0.95)
lines(COD.x.pred.R$Nsc, COD.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(COD.x.pred.R$Nsc, COD.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
COD.x.pred.G <- data.frame(COD.lNsc = seq(min(COD.lNsc,na.rm=T), max(COD.lNsc,na.rm=T), length.out=100))
COD.pred.interval.G <- predict(COD.gomp, newdata=COD.x.pred.G, interval="confidence", level=0.95)
lines(exp(COD.x.pred.G$COD.lNsc), COD.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(COD.x.pred.G$COD.lNsc), COD.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
COD.R.N.descaled <- attr(COD.Nsc, 'scaled:scale') * COD.x.pred.R$Nsc
COD.G.N.descaled <- attr(COD.Nsc, 'scaled:scale') * exp(COD.x.pred.G$COD.lNsc)

## export results
COD.results.out <- data.frame(x.R=COD.R.N.descaled,
                              R.fit=COD.pred.interval.R[,"fit"],
                              R.upr=COD.pred.interval.R[,"upr"],
                              R.lwr=COD.pred.interval.R[,"lwr"],
                              x.G=COD.G.N.descaled,
                              G.fit=COD.pred.interval.G[,"fit"],
                              G.upr=COD.pred.interval.G[,"upr"],
                              G.lwr=COD.pred.interval.G[,"lwr"])
head(COD.results.out)

# save to .csv
write.csv(COD.results.out, "COD_Ricker_Gompertz_fits.csv", row.names=F)


# Mali
MLI <- read.csv("MLI.csv", header=T)
head(MLI)
MLI.Nsc <- scale(MLI[,2], scale=T, center=F)
MLI$Nsc <- as.numeric(MLI.Nsc)
MLI$r <- c(log(MLI$Nsc[2:dim(popreg)[1]] / MLI$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(MLI)
plot(MLI$Nsc, MLI$r, xlab="N", ylab="r", pch=19)
MLI.use <- MLI[which(MLI$year >= 2004),]
plot(MLI.use$Nsc, MLI.use$r, xlab="N", ylab="r", pch=19)

# Ricker
MLI.rick <- lm(r ~ Nsc, data=MLI.use)
abline(MLI.rick, lty=2, col="red")
summary(MLI.rick)

MLI.acf.pacf <- plot_combined_acf_pacf(na.omit(MLI$r), max_lag = 40)
MLI.acf.pacf

## Newey-West standard errors
MLI.coeftest <- coeftest(MLI.rick, vcov = NeweyWest(MLI.rick, prewhite = FALSE), save=T)
MLI.coeftest
MLI.coeftest.confint <- confint(MLI.coeftest)
MLI.coeftest.confint

MLI.coeftest.fitted <- attr(MLI.coeftest, "object")[[5]]
MLI.R2.NW <- cor(na.omit(MLI.use$r), MLI.coeftest.fitted)^2
MLI.R2.NW

# Gompertz
MLI.lNsc <- log(MLI.use$Nsc)
MLI.gomp <- lm(MLI.use$r ~ MLI.lNsc)
summary(MLI.gomp)

MLI.RLGL2.AIC.vec <- c(AICc(MLI.rick), AICc(MLI.gomp))
MLI.RLGL2.dAIC.vec <- delta.IC(MLI.RLGL2.AIC.vec)
MLI.RLGL2.wAIC.vec <- weight.IC(MLI.RLGL2.dAIC.vec)
print(MLI.RLGL2.wAIC.vec)
MLI.ER.RLGL2 <- MLI.RLGL2.wAIC.vec[1]/MLI.RLGL2.wAIC.vec[2]
MLI.ER.RLGL2
1/MLI.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(MLI.use$Nsc, MLI.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
MLI.x.pred.R <- data.frame(Nsc = seq(min(MLI.use$Nsc,na.rm=T), max(MLI.use$Nsc,na.rm=T), length.out=100))
MLI.pred.interval.R <- predict(MLI.rick, newdata=MLI.x.pred.R, interval="confidence", level=0.95)
lines(MLI.x.pred.R$Nsc, MLI.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(MLI.x.pred.R$Nsc, MLI.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
MLI.x.pred.G <- data.frame(MLI.lNsc = seq(min(MLI.lNsc,na.rm=T), max(MLI.lNsc,na.rm=T), length.out=100))
MLI.pred.interval.G <- predict(MLI.gomp, newdata=MLI.x.pred.G, interval="confidence", level=0.95)
lines(exp(MLI.x.pred.G$MLI.lNsc), MLI.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(MLI.x.pred.G$MLI.lNsc), MLI.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
MLI.R.N.descaled <- attr(MLI.Nsc, 'scaled:scale') * MLI.x.pred.R$Nsc
MLI.G.N.descaled <- attr(MLI.Nsc, 'scaled:scale') * exp(MLI.x.pred.G$MLI.lNsc)

## export results
MLI.results.out <- data.frame(x.R=MLI.R.N.descaled,
                              R.fit=MLI.pred.interval.R[,"fit"],
                              R.upr=MLI.pred.interval.R[,"upr"],
                              R.lwr=MLI.pred.interval.R[,"lwr"],
                              x.G=MLI.G.N.descaled,
                              G.fit=MLI.pred.interval.G[,"fit"],
                              G.upr=MLI.pred.interval.G[,"upr"],
                              G.lwr=MLI.pred.interval.G[,"lwr"])
head(MLI.results.out)

# save to .csv
write.csv(MLI.results.out, "MLI_Ricker_Gompertz_fits.csv", row.names=F)
                                                        
                                                        
# Chad
TCD <- read.csv("TCD.csv", header=T)
head(TCD)
TCD.Nsc <- scale(TCD[,2], scale=T, center=F)
TCD$Nsc <- as.numeric(TCD.Nsc)
TCD$r <- c(log(TCD$Nsc[2:dim(popreg)[1]] / TCD$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(TCD)
plot(TCD$Nsc, TCD$r, xlab="N", ylab="r", pch=19)
TCD.use <- TCD[which(TCD$year >= 2003),]
plot(TCD.use$Nsc, TCD.use$r, xlab="N", ylab="r", pch=19)

# Ricker
TCD.rick <- lm(r ~ Nsc, data=TCD.use)
abline(TCD.rick, lty=2, col="red")
summary(TCD.rick)

TCD.acf.pacf <- plot_combined_acf_pacf(na.omit(TCD$r), max_lag = 40)
TCD.acf.pacf

## Newey-West standard errors
TCD.coeftest <- coeftest(TCD.rick, vcov = NeweyWest(TCD.rick, prewhite = FALSE), save=T)
TCD.coeftest
TCD.coeftest.confint <- confint(TCD.coeftest)
TCD.coeftest.confint

TCD.coeftest.fitted <- attr(TCD.coeftest, "object")[[5]]
TCD.R2.NW <- cor(na.omit(TCD.use$r), TCD.coeftest.fitted)^2
TCD.R2.NW

# Gompertz
TCD.lNsc <- log(TCD.use$Nsc)
TCD.gomp <- lm(TCD.use$r ~ TCD.lNsc)
summary(TCD.gomp)

TCD.RLGL2.AIC.vec <- c(AICc(TCD.rick), AICc(TCD.gomp))
TCD.RLGL2.dAIC.vec <- delta.IC(TCD.RLGL2.AIC.vec)
TCD.RLGL2.wAIC.vec <- weight.IC(TCD.RLGL2.dAIC.vec)
print(TCD.RLGL2.wAIC.vec)
TCD.ER.RLGL2 <- TCD.RLGL2.wAIC.vec[1]/TCD.RLGL2.wAIC.vec[2]
TCD.ER.RLGL2
1/TCD.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(TCD.use$Nsc, TCD.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
TCD.x.pred.R <- data.frame(Nsc = seq(min(TCD.use$Nsc,na.rm=T), max(TCD.use$Nsc,na.rm=T), length.out=100))
TCD.pred.interval.R <- predict(TCD.rick, newdata=TCD.x.pred.R, interval="confidence", level=0.95)
lines(TCD.x.pred.R$Nsc, TCD.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(TCD.x.pred.R$Nsc, TCD.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
TCD.x.pred.G <- data.frame(TCD.lNsc = seq(min(TCD.lNsc,na.rm=T), max(TCD.lNsc,na.rm=T), length.out=100))
TCD.pred.interval.G <- predict(TCD.gomp, newdata=TCD.x.pred.G, interval="confidence", level=0.95)
lines(exp(TCD.x.pred.G$TCD.lNsc), TCD.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(TCD.x.pred.G$TCD.lNsc), TCD.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
TCD.R.N.descaled <- attr(TCD.Nsc, 'scaled:scale') * TCD.x.pred.R$Nsc
TCD.G.N.descaled <- attr(TCD.Nsc, 'scaled:scale') * exp(TCD.x.pred.G$TCD.lNsc)

## export results
TCD.results.out <- data.frame(x.R=TCD.R.N.descaled,
                              R.fit=TCD.pred.interval.R[,"fit"],
                              R.upr=TCD.pred.interval.R[,"upr"],
                              R.lwr=TCD.pred.interval.R[,"lwr"],
                              x.G=TCD.G.N.descaled,
                              G.fit=TCD.pred.interval.G[,"fit"],
                              G.upr=TCD.pred.interval.G[,"upr"],
                              G.lwr=TCD.pred.interval.G[,"lwr"])
head(TCD.results.out)

# save to .csv
write.csv(TCD.results.out, "TCD_Ricker_Gompertz_fits.csv", row.names=F)


# Angola
AGO <- read.csv("AGO.csv", header=T)
head(AGO)
AGO.Nsc <- scale(AGO[,2], scale=T, center=F)
AGO$Nsc <- as.numeric(AGO.Nsc)
AGO$r <- c(log(AGO$Nsc[2:dim(popreg)[1]] / AGO$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(AGO)
plot(AGO$Nsc, AGO$r, xlab="N", ylab="r", pch=19)
AGO.use <- AGO[which(AGO$year >= 2010),]
plot(AGO.use$Nsc, AGO.use$r, xlab="N", ylab="r", pch=19)

# Ricker
AGO.rick <- lm(r ~ Nsc, data=AGO.use)
abline(AGO.rick, lty=2, col="red")
summary(AGO.rick)

AGO.acf.pacf <- plot_combined_acf_pacf(na.omit(AGO$r), max_lag = 40)
AGO.acf.pacf

## Newey-West standard errors
AGO.coeftest <- coeftest(AGO.rick, vcov = NeweyWest(AGO.rick, prewhite = FALSE), save=T)
AGO.coeftest
AGO.coeftest.confint <- confint(AGO.coeftest)
AGO.coeftest.confint

AGO.coeftest.fitted <- attr(AGO.coeftest, "object")[[5]]
AGO.R2.NW <- cor(na.omit(AGO.use$r), AGO.coeftest.fitted)^2
AGO.R2.NW

# Gompertz
AGO.lNsc <- log(AGO.use$Nsc)
AGO.gomp <- lm(AGO.use$r ~ AGO.lNsc)
summary(AGO.gomp)

AGO.RLGL2.AIC.vec <- c(AICc(AGO.rick), AICc(AGO.gomp))
AGO.RLGL2.dAIC.vec <- delta.IC(AGO.RLGL2.AIC.vec)
AGO.RLGL2.wAIC.vec <- weight.IC(AGO.RLGL2.dAIC.vec)
print(AGO.RLGL2.wAIC.vec)
AGO.ER.RLGL2 <- AGO.RLGL2.wAIC.vec[1]/AGO.RLGL2.wAIC.vec[2]
AGO.ER.RLGL2
1/AGO.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(AGO.use$Nsc, AGO.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
AGO.x.pred.R <- data.frame(Nsc = seq(min(AGO.use$Nsc,na.rm=T), max(AGO.use$Nsc,na.rm=T), length.out=100))
AGO.pred.interval.R <- predict(AGO.rick, newdata=AGO.x.pred.R, interval="confidence", level=0.95)
lines(AGO.x.pred.R$Nsc, AGO.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(AGO.x.pred.R$Nsc, AGO.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
AGO.x.pred.G <- data.frame(AGO.lNsc = seq(min(AGO.lNsc,na.rm=T), max(AGO.lNsc,na.rm=T), length.out=100))
AGO.pred.interval.G <- predict(AGO.gomp, newdata=AGO.x.pred.G, interval="confidence", level=0.95)
lines(exp(AGO.x.pred.G$AGO.lNsc), AGO.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(AGO.x.pred.G$AGO.lNsc), AGO.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
AGO.R.N.descaled <- attr(AGO.Nsc, 'scaled:scale') * AGO.x.pred.R$Nsc
AGO.G.N.descaled <- attr(AGO.Nsc, 'scaled:scale') * exp(AGO.x.pred.G$AGO.lNsc)  

## export results
AGO.results.out <- data.frame(x.R=AGO.R.N.descaled,
                              R.fit=AGO.pred.interval.R[,"fit"],
                              R.upr=AGO.pred.interval.R[,"upr"],
                              R.lwr=AGO.pred.interval.R[,"lwr"],
                              x.G=AGO.G.N.descaled,
                              G.fit=AGO.pred.interval.G[,"fit"],
                              G.upr=AGO.pred.interval.G[,"upr"],
                              G.lwr=AGO.pred.interval.G[,"lwr"])
head(AGO.results.out)

# save to .csv
write.csv(AGO.results.out, "AGO_Ricker_Gompertz_fits.csv", row.names=F)


# Nigeria
NGA <- read.csv("NGA.csv", header=T)
head(NGA)
NGA.Nsc <- scale(NGA[,2], scale=T, center=F)
NGA$Nsc <- as.numeric(NGA.Nsc)
NGA$r <- c(log(NGA$Nsc[2:dim(popreg)[1]] / NGA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(NGA)
plot(NGA$Nsc, NGA$r, xlab="N", ylab="r", pch=19)
NGA.use <- NGA[which(NGA$year >= 2010),]
plot(NGA.use$Nsc, NGA.use$r, xlab="N", ylab="r", pch=19)

# Ricker
NGA.rick <- lm(r ~ Nsc, data=NGA.use)
abline(NGA.rick, lty=2, col="red")
summary(NGA.rick)

NGA.acf.pacf <- plot_combined_acf_pacf(na.omit(NGA$r), max_lag = 40)
NGA.acf.pacf

## Newey-West standard errors
NGA.coeftest <- coeftest(NGA.rick, vcov = NeweyWest(NGA.rick, prewhite = FALSE), save=T)
NGA.coeftest
NGA.coeftest.confint <- confint(NGA.coeftest)
NGA.coeftest.confint

NGA.coeftest.fitted <- attr(NGA.coeftest, "object")[[5]]
NGA.R2.NW <- cor(na.omit(NGA.use$r), NGA.coeftest.fitted)^2
NGA.R2.NW

# Gompertz
NGA.lNsc <- log(NGA.use$Nsc)
NGA.gomp <- lm(NGA.use$r ~ NGA.lNsc)
summary(NGA.gomp)

NGA.RLGL2.AIC.vec <- c(AICc(NGA.rick), AICc(NGA.gomp))
NGA.RLGL2.dAIC.vec <- delta.IC(NGA.RLGL2.AIC.vec)
NGA.RLGL2.wAIC.vec <- weight.IC(NGA.RLGL2.dAIC.vec)
print(NGA.RLGL2.wAIC.vec)
NGA.ER.RLGL2 <- NGA.RLGL2.wAIC.vec[1]/NGA.RLGL2.wAIC.vec[2]
NGA.ER.RLGL2
1/NGA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(NGA.use$Nsc, NGA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
NGA.x.pred.R <- data.frame(Nsc = seq(min(NGA.use$Nsc,na.rm=T), max(NGA.use$Nsc,na.rm=T), length.out=100))
NGA.pred.interval.R <- predict(NGA.rick, newdata=NGA.x.pred.R, interval="confidence", level=0.95)
lines(NGA.x.pred.R$Nsc, NGA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(NGA.x.pred.R$Nsc, NGA.pred.interval.R[,"upr"], col="blue", lty=2)
     
# Gompertz
NGA.x.pred.G <- data.frame(NGA.lNsc = seq(min(NGA.lNsc,na.rm=T), max(NGA.lNsc,na.rm=T), length.out=100))
NGA.pred.interval.G <- predict(NGA.gomp, newdata=NGA.x.pred.G, interval="confidence", level=0.95)
lines(exp(NGA.x.pred.G$NGA.lNsc), NGA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(NGA.x.pred.G$NGA.lNsc), NGA.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
NGA.R.N.descaled <- attr(NGA.Nsc, 'scaled:scale') * NGA.x.pred.R$Nsc
NGA.G.N.descaled <- attr(NGA.Nsc, 'scaled:scale') * exp(NGA.x.pred.G$NGA.lNsc)

## export results
NGA.results.out <- data.frame(x.R=NGA.R.N.descaled,
                              R.fit=NGA.pred.interval.R[,"fit"],
                              R.upr=NGA.pred.interval.R[,"upr"],
                              R.lwr=NGA.pred.interval.R[,"lwr"],
                              x.G=NGA.G.N.descaled,
                              G.fit=NGA.pred.interval.G[,"fit"],
                              G.upr=NGA.pred.interval.G[,"upr"],
                              G.lwr=NGA.pred.interval.G[,"lwr"])
head(NGA.results.out)

# save to .csv
write.csv(NGA.results.out, "NGA_Ricker_Gompertz_fits.csv", row.names=F)


# Burundi
BDI <- read.csv("BDI.csv", header=T)
head(BDI)
BDI.Nsc <- scale(BDI[,2], scale=T, center=F)
BDI$Nsc <- as.numeric(BDI.Nsc)
BDI$r <- c(log(BDI$Nsc[2:dim(popreg)[1]] / BDI$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(BDI)
plot(BDI$Nsc, BDI$r, xlab="N", ylab="r", pch=19)
BDI.use <- BDI[which(BDI$year >= 2009),]
plot(BDI.use$Nsc, BDI.use$r, xlab="N", ylab="r", pch=19)

# Ricker
BDI.rick <- lm(r ~ Nsc, data=BDI.use)
abline(BDI.rick, lty=2, col="red")
summary(BDI.rick)

BDI.acf.pacf <- plot_combined_acf_pacf(na.omit(BDI$r), max_lag = 40)
BDI.acf.pacf

## Newey-West standard errors
BDI.coeftest <- coeftest(BDI.rick, vcov = NeweyWest(BDI.rick, prewhite = FALSE), save=T)
BDI.coeftest
BDI.coeftest.confint <- confint(BDI.coeftest)
BDI.coeftest.confint

BDI.coeftest.fitted <- attr(BDI.coeftest, "object")[[5]]
BDI.R2.NW <- cor(na.omit(BDI.use$r), BDI.coeftest.fitted)^2
BDI.R2.NW

# Gompertz
BDI.lNsc <- log(BDI.use$Nsc)
BDI.gomp <- lm(BDI.use$r ~ BDI.lNsc)
summary(BDI.gomp)

BDI.RLGL2.AIC.vec <- c(AICc(BDI.rick), AICc(BDI.gomp))
BDI.RLGL2.dAIC.vec <- delta.IC(BDI.RLGL2.AIC.vec)
BDI.RLGL2.wAIC.vec <- weight.IC(BDI.RLGL2.dAIC.vec)
print(BDI.RLGL2.wAIC.vec)
BDI.ER.RLGL2 <- BDI.RLGL2.wAIC.vec[1]/BDI.RLGL2.wAIC.vec[2]
BDI.ER.RLGL2
1/BDI.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(BDI.use$Nsc, BDI.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
BDI.x.pred.R <- data.frame(Nsc = seq(min(BDI.use$Nsc,na.rm=T), max(BDI.use$Nsc,na.rm=T), length.out=100))
BDI.pred.interval.R <- predict(BDI.rick, newdata=BDI.x.pred.R, interval="confidence", level=0.95)
lines(BDI.x.pred.R$Nsc, BDI.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(BDI.x.pred.R$Nsc, BDI.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
BDI.x.pred.G <- data.frame(BDI.lNsc = seq(min(BDI.lNsc,na.rm=T), max(BDI.lNsc,na.rm=T), length.out=100))
BDI.pred.interval.G <- predict(BDI.gomp, newdata=BDI.x.pred.G, interval="confidence", level=0.95)
lines(exp(BDI.x.pred.G$BDI.lNsc), BDI.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(BDI.x.pred.G$BDI.lNsc), BDI.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
BDI.R.N.descaled <- attr(BDI.Nsc, 'scaled:scale') * BDI.x.pred.R$Nsc
BDI.G.N.descaled <- attr(BDI.Nsc, 'scaled:scale') * exp(BDI.x.pred.G$BDI.lNsc)

## export results
BDI.results.out <- data.frame(x.R=BDI.R.N.descaled,
                              R.fit=BDI.pred.interval.R[,"fit"],
                              R.upr=BDI.pred.interval.R[,"upr"],
                              R.lwr=BDI.pred.interval.R[,"lwr"],
                              x.G=BDI.G.N.descaled,
                              G.fit=BDI.pred.interval.G[,"fit"],
                              G.upr=BDI.pred.interval.G[,"upr"],
                              G.lwr=BDI.pred.interval.G[,"lwr"])
head(BDI.results.out)

# save to .csv
write.csv(BDI.results.out, "BDI_Ricker_Gompertz_fits.csv", row.names=F)


# Burkina Faso
BFA <- read.csv("BFA.csv", header=T)
head(BFA)
BFA.Nsc <- scale(BFA[,2], scale=T, center=F)
BFA$Nsc <- as.numeric(BFA.Nsc)
BFA$r <- c(log(BFA$Nsc[2:dim(popreg)[1]] / BFA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(BFA)
plot(BFA$Nsc, BFA$r, xlab="N", ylab="r", pch=19)
BFA.use <- BFA[which(BFA$year >= 2004),]
plot(BFA.use$Nsc, BFA.use$r, xlab="N", ylab="r", pch=19)

# Ricker
BFA.rick <- lm(r ~ Nsc, data=BFA.use)
abline(BFA.rick, lty=2, col="red")
summary(BFA.rick)

BFA.acf.pacf <- plot_combined_acf_pacf(na.omit(BFA$r), max_lag = 40)
BFA.acf.pacf

## Newey-West standard errors
BFA.coeftest <- coeftest(BFA.rick, vcov = NeweyWest(BFA.rick, prewhite = FALSE), save=T)
BFA.coeftest
BFA.coeftest.confint <- confint(BFA.coeftest)
BFA.coeftest.confint

BFA.coeftest.fitted <- attr(BFA.coeftest, "object")[[5]]
BFA.R2.NW <- cor(na.omit(BFA.use$r), BFA.coeftest.fitted)^2
BFA.R2.NW

# Gompertz
BFA.lNsc <- log(BFA.use$Nsc)
BFA.gomp <- lm(BFA.use$r ~ BFA.lNsc)
summary(BFA.gomp)

BFA.RLGL2.AIC.vec <- c(AICc(BFA.rick), AICc(BFA.gomp))
BFA.RLGL2.dAIC.vec <- delta.IC(BFA.RLGL2.AIC.vec)
BFA.RLGL2.wAIC.vec <- weight.IC(BFA.RLGL2.dAIC.vec)
print(BFA.RLGL2.wAIC.vec)
BFA.ER.RLGL2 <- BFA.RLGL2.wAIC.vec[1]/BFA.RLGL2.wAIC.vec[2]
BFA.ER.RLGL2
1/BFA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(BFA.use$Nsc, BFA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
BFA.x.pred.R <- data.frame(Nsc = seq(min(BFA.use$Nsc,na.rm=T), max(BFA.use$Nsc,na.rm=T), length.out=100))
BFA.pred.interval.R <- predict(BFA.rick, newdata=BFA.x.pred.R, interval="confidence", level=0.95)
lines(BFA.x.pred.R$Nsc, BFA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(BFA.x.pred.R$Nsc, BFA.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
BFA.x.pred.G <- data.frame(BFA.lNsc = seq(min(BFA.lNsc,na.rm=T), max(BFA.lNsc,na.rm=T), length.out=100))
BFA.pred.interval.G <- predict(BFA.gomp, newdata=BFA.x.pred.G, interval="confidence", level=0.95)
lines(exp(BFA.x.pred.G$BFA.lNsc), BFA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(BFA.x.pred.G$BFA.lNsc), BFA.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
BFA.R.N.descaled <- attr(BFA.Nsc, 'scaled:scale') * BFA.x.pred.R$Nsc
BFA.G.N.descaled <- attr(BFA.Nsc, 'scaled:scale') * exp(BFA.x.pred.G$BFA.lNsc)

## export results
BFA.results.out <- data.frame(x.R=BFA.R.N.descaled,
                              R.fit=BFA.pred.interval.R[,"fit"],
                              R.upr=BFA.pred.interval.R[,"upr"],
                              R.lwr=BFA.pred.interval.R[,"lwr"],
                              x.G=BFA.G.N.descaled,
                              G.fit=BFA.pred.interval.G[,"fit"],
                              G.upr=BFA.pred.interval.G[,"upr"],
                              G.lwr=BFA.pred.interval.G[,"lwr"])
head(BFA.results.out)

# save to .csv
write.csv(BFA.results.out, "BFA_Ricker_Gompertz_fits.csv", row.names=F)



# Gambia
GMB <- read.csv("GMB.csv", header=T)
head(GMB)
GMB.Nsc <- scale(GMB[,2], scale=T, center=F)
GMB$Nsc <- as.numeric(GMB.Nsc)
GMB$r <- c(log(GMB$Nsc[2:dim(popreg)[1]] / GMB$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(GMB)
plot(GMB$Nsc, GMB$r, xlab="N", ylab="r", pch=19)
GMB.use <- GMB[which(GMB$year >= 2012),]
plot(GMB.use$Nsc, GMB.use$r, xlab="N", ylab="r", pch=19)

# Ricker
GMB.rick <- lm(r ~ Nsc, data=GMB.use)
abline(GMB.rick, lty=2, col="red")
summary(GMB.rick)

GMB.acf.pacf <- plot_combined_acf_pacf(na.omit(GMB$r), max_lag = 40)
GMB.acf.pacf

## Newey-West standard errors
GMB.coeftest <- coeftest(GMB.rick, vcov = NeweyWest(GMB.rick, prewhite = FALSE), save=T)
GMB.coeftest
GMB.coeftest.confint <- confint(GMB.coeftest)
GMB.coeftest.confint

GMB.coeftest.fitted <- attr(GMB.coeftest, "object")[[5]]
GMB.R2.NW <- cor(na.omit(GMB.use$r), GMB.coeftest.fitted)^2
GMB.R2.NW

# Gompertz
GMB.lNsc <- log(GMB.use$Nsc)
GMB.gomp <- lm(GMB.use$r ~ GMB.lNsc)
summary(GMB.gomp)

GMB.RLGL2.AIC.vec <- c(AICc(GMB.rick), AICc(GMB.gomp))
GMB.RLGL2.dAIC.vec <- delta.IC(GMB.RLGL2.AIC.vec)
GMB.RLGL2.wAIC.vec <- weight.IC(GMB.RLGL2.dAIC.vec)
print(GMB.RLGL2.wAIC.vec)
GMB.ER.RLGL2 <- GMB.RLGL2.wAIC.vec[1]/GMB.RLGL2.wAIC.vec[2]
GMB.ER.RLGL2
1/GMB.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(GMB.use$Nsc, GMB.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
GMB.x.pred.R <- data.frame(Nsc = seq(min(GMB.use$Nsc,na.rm=T), max(GMB.use$Nsc,na.rm=T), length.out=100))
GMB.pred.interval.R <- predict(GMB.rick, newdata=GMB.x.pred.R, interval="confidence", level=0.95)
lines(GMB.x.pred.R$Nsc, GMB.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(GMB.x.pred.R$Nsc, GMB.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
GMB.x.pred.G <- data.frame(GMB.lNsc = seq(min(GMB.lNsc,na.rm=T), max(GMB.lNsc,na.rm=T), length.out=100))
GMB.pred.interval.G <- predict(GMB.gomp, newdata=GMB.x.pred.G, interval="confidence", level=0.95)
lines(exp(GMB.x.pred.G$GMB.lNsc), GMB.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(GMB.x.pred.G$GMB.lNsc), GMB.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
GMB.R.N.descaled <- attr(GMB.Nsc, 'scaled:scale') * GMB.x.pred.R$Nsc
GMB.G.N.descaled <- attr(GMB.Nsc, 'scaled:scale') * exp(GMB.x.pred.G$GMB.lNsc)

## export results
GMB.results.out <- data.frame(x.R=GMB.R.N.descaled,
                              R.fit=GMB.pred.interval.R[,"fit"],
                              R.upr=GMB.pred.interval.R[,"upr"],
                              R.lwr=GMB.pred.interval.R[,"lwr"],
                              x.G=GMB.G.N.descaled,
                              G.fit=GMB.pred.interval.G[,"fit"],
                              G.upr=GMB.pred.interval.G[,"upr"],
                              G.lwr=GMB.pred.interval.G[,"lwr"])
head(GMB.results.out)

# save to .csv
write.csv(GMB.results.out, "GMB_Ricker_Gompertz_fits.csv", row.names=F)



# Uganda
UGA <- read.csv("UGA.csv", header=T)
head(UGA)
UGA.Nsc <- scale(UGA[,2], scale=T, center=F)
UGA$Nsc <- as.numeric(UGA.Nsc)
UGA$r <- c(log(UGA$Nsc[2:dim(popreg)[1]] / UGA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(UGA)
plot(UGA$Nsc, UGA$r, xlab="N", ylab="r", pch=19)
UGA.use <- UGA[which(UGA$year >= 2016),]
plot(UGA.use$Nsc, UGA.use$r, xlab="N", ylab="r", pch=19)

# Ricker
UGA.rick <- lm(r ~ Nsc, data=UGA.use)
abline(UGA.rick, lty=2, col="red")
summary(UGA.rick)

UGA.acf.pacf <- plot_combined_acf_pacf(na.omit(UGA$r), max_lag = 40)
UGA.acf.pacf

## Newey-West standard errors
UGA.coeftest <- coeftest(UGA.rick, vcov = NeweyWest(UGA.rick, prewhite = FALSE), save=T)
UGA.coeftest
UGA.coeftest.confint <- confint(UGA.coeftest)
UGA.coeftest.confint

UGA.coeftest.fitted <- attr(UGA.coeftest, "object")[[5]]
UGA.R2.NW <- cor(na.omit(UGA.use$r), UGA.coeftest.fitted)^2
UGA.R2.NW

# Gompertz
UGA.lNsc <- log(UGA.use$Nsc)
UGA.gomp <- lm(UGA.use$r ~ UGA.lNsc)
summary(UGA.gomp)

UGA.RLGL2.AIC.vec <- c(AICc(UGA.rick), AICc(UGA.gomp))
UGA.RLGL2.dAIC.vec <- delta.IC(UGA.RLGL2.AIC.vec)
UGA.RLGL2.wAIC.vec <- weight.IC(UGA.RLGL2.dAIC.vec)
print(UGA.RLGL2.wAIC.vec)
UGA.ER.RLGL2 <- UGA.RLGL2.wAIC.vec[1]/UGA.RLGL2.wAIC.vec[2]
UGA.ER.RLGL2
1/UGA.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(UGA.use$Nsc, UGA.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
UGA.x.pred.R <- data.frame(Nsc = seq(min(UGA.use$Nsc,na.rm=T), max(UGA.use$Nsc,na.rm=T), length.out=100))
UGA.pred.interval.R <- predict(UGA.rick, newdata=UGA.x.pred.R, interval="confidence", level=0.95)
lines(UGA.x.pred.R$Nsc, UGA.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(UGA.x.pred.R$Nsc, UGA.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
UGA.x.pred.G <- data.frame(UGA.lNsc = seq(min(UGA.lNsc,na.rm=T), max(UGA.lNsc,na.rm=T), length.out=100))
UGA.pred.interval.G <- predict(UGA.gomp, newdata=UGA.x.pred.G, interval="confidence", level=0.95)
lines(exp(UGA.x.pred.G$UGA.lNsc), UGA.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(UGA.x.pred.G$UGA.lNsc), UGA.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
UGA.R.N.descaled <- attr(UGA.Nsc, 'scaled:scale') * UGA.x.pred.R$Nsc
UGA.G.N.descaled <- attr(UGA.Nsc, 'scaled:scale') * exp(UGA.x.pred.G$UGA.lNsc)

## export results
UGA.results.out <- data.frame(x.R=UGA.R.N.descaled,
                              R.fit=UGA.pred.interval.R[,"fit"],
                              R.upr=UGA.pred.interval.R[,"upr"],
                              R.lwr=UGA.pred.interval.R[,"lwr"],
                              x.G=UGA.G.N.descaled,
                              G.fit=UGA.pred.interval.G[,"fit"],
                              G.upr=UGA.pred.interval.G[,"upr"],
                              G.lwr=UGA.pred.interval.G[,"lwr"])
head(UGA.results.out)

# save to .csv
write.csv(UGA.results.out, "UGA_Ricker_Gompertz_fits.csv", row.names=F)


## China
CHN <- read.csv("china1950-2023.csv", header=T)
head(CHN)
CHN.Nsc <- scale(CHN[,2], scale=T, center=F)
CHN$Nsc <- as.numeric(CHN.Nsc)
CHN$r <- c(log(CHN$Nsc[2:dim(popreg)[1]] / CHN$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(CHN)
plot(CHN$Nsc, CHN$r, xlab="N", ylab="r", pch=19)
CHN.use <- CHN[which(CHN$year >= 1962),]
plot(CHN.use$Nsc, CHN.use$r, xlab="N", ylab="r", pch=19)

# Ricker
CHN.rick <- lm(r ~ Nsc, data=CHN.use)
abline(CHN.rick, lty=2, col="red")
summary(CHN.rick)

CHN.acf.pacf <- plot_combined_acf_pacf(na.omit(CHN$r), max_lag = 40)
CHN.acf.pacf

## Newey-West standard errors
CHN.coeftest <- coeftest(CHN.rick, vcov = NeweyWest(CHN.rick, prewhite = FALSE), save=T)
CHN.coeftest
CHN.coeftest.confint <- confint(CHN.coeftest)
CHN.coeftest.confint

CHN.coeftest.fitted <- attr(CHN.coeftest, "object")[[5]]
CHN.R2.NW <- cor(na.omit(CHN.use$r), CHN.coeftest.fitted)^2
CHN.R2.NW

# Gompertz
CHN.lNsc <- log(CHN.use$Nsc)
CHN.gomp <- lm(CHN.use$r ~ CHN.lNsc)
summary(CHN.gomp)

CHN.RLGL2.AIC.vec <- c(AICc(CHN.rick), AICc(CHN.gomp))
CHN.RLGL2.dAIC.vec <- delta.IC(CHN.RLGL2.AIC.vec)
CHN.RLGL2.wAIC.vec <- weight.IC(CHN.RLGL2.dAIC.vec)
print(CHN.RLGL2.wAIC.vec)
CHN.ER.RLGL2 <- CHN.RLGL2.wAIC.vec[1]/CHN.RLGL2.wAIC.vec[2]
CHN.ER.RLGL2
1/CHN.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(CHN.use$Nsc, CHN.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
CHN.x.pred.R <- data.frame(Nsc = seq(min(CHN.use$Nsc,na.rm=T), max(CHN.use$Nsc,na.rm=T), length.out=100))
CHN.pred.interval.R <- predict(CHN.rick, newdata=CHN.x.pred.R, interval="confidence", level=0.95)
lines(CHN.x.pred.R$Nsc, CHN.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(CHN.x.pred.R$Nsc, CHN.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
CHN.x.pred.G <- data.frame(CHN.lNsc = seq(min(CHN.lNsc,na.rm=T), max(CHN.lNsc,na.rm=T), length.out=100))
CHN.pred.interval.G <- predict(CHN.gomp, newdata=CHN.x.pred.G, interval="confidence", level=0.95)
lines(exp(CHN.x.pred.G$CHN.lNsc), CHN.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(CHN.x.pred.G$CHN.lNsc), CHN.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
CHN.R.N.descaled <- attr(CHN.Nsc, 'scaled:scale') * CHN.x.pred.R$Nsc
CHN.G.N.descaled <- attr(CHN.Nsc, 'scaled:scale') * exp(CHN.x.pred.G$CHN.lNsc)

## export results
CHN.results.out <- data.frame(x.R=CHN.R.N.descaled,
                              R.fit=CHN.pred.interval.R[,"fit"],
                              R.upr=CHN.pred.interval.R[,"upr"],
                              R.lwr=CHN.pred.interval.R[,"lwr"],
                              x.G=CHN.G.N.descaled,
                              G.fit=CHN.pred.interval.G[,"fit"],
                              G.upr=CHN.pred.interval.G[,"upr"],
                              G.lwr=CHN.pred.interval.G[,"lwr"])
head(CHN.results.out)

# save to .csv
write.csv(CHN.results.out, "CHN_Ricker_Gompertz_fits.csv", row.names=F)




# EAST AND SOUTH-EASTERN ASIA EXCLUDING CHINA
ESEA.nCHN <- ESEA
head(ESEA.nCHN)
head(CHN)

ESEA.nCHN$ESEA <- ESEA.nCHN$ESEA - CHN$N # remove China
head(ESEA.nCHN)

ESEA.nCHN.NSc <- scale(ESEA.nCHN[,2], scale=T, center=F)
ESEA.nCHN$Nsc <- as.numeric(ESEA.nCHN.NSc)
ESEA.nCHN$r <- c(log(ESEA.nCHN$Nsc[2:dim(popreg)[1]] / ESEA.nCHN$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(ESEA.nCHN)
head(ESEA)
plot(ESEA.nCHN$Nsc, ESEA.nCHN$r, xlab="N", ylab="r", pch=19, type="b")

# export
write.csv(ESEA.nCHN, "ESEA_minus_CHN.csv", row.names=F)

year.neg.st <- ESEA.nCHN$year[which(ESEA.nCHN$r == max(ESEA.nCHN$r, na.rm=T))]
ESEA.nCHN.use <- ESEA.nCHN[which(ESEA.nCHN$year >= year.neg.st),]
plot(ESEA.nCHN.use$Nsc, ESEA.nCHN.use$r, xlab="N", ylab="r", pch=19)

# plot all ESEA for comparison
plot(ESEA.use$Nsc, ESEA.use$r, xlab="N", ylab="r", pch=19)
plot(ESEA$ESEA, ESEA.nCHN$ESEA, xlab="ESEA total", ylab="ESEA minus China", pch=19)

# Ricker
plot(ESEA.nCHN.use$Nsc, ESEA.nCHN.use$r, xlab="N", ylab="r", pch=19)
ESEA.nCHN.rick <- lm(r ~ Nsc, data=ESEA.nCHN.use)
abline(ESEA.nCHN.rick, lty=2, col="red")
summary(ESEA.nCHN.rick)

ESEA.nCHN.acf.pacf <- plot_combined_acf_pacf(na.omit(ESEA.nCHN$r), max_lag = 40)
ESEA.nCHN.acf.pacf

## Newey-West standard errors
ESEA.nCHN.coeftest <- coeftest(ESEA.nCHN.rick, vcov = NeweyWest(ESEA.nCHN.rick, prewhite = FALSE), save=T)
ESEA.nCHN.coeftest
ESEA.nCHN.coeftest.confint <- confint(ESEA.nCHN.coeftest)
ESEA.nCHN.coeftest.confint

ESEA.nCHN.coeftest.fitted <- attr(ESEA.nCHN.coeftest, "object")[[5]]
ESEA.nCHN.R2.NW <- cor(na.omit(ESEA.nCHN.use$r), ESEA.nCHN.coeftest.fitted)^2
ESEA.nCHN.R2.NW

# Gompertz
ESEA.nCHN.lNsc <- log(ESEA.nCHN.use$Nsc)
ESEA.nCHN.gomp <- lm(ESEA.nCHN.use$r ~ ESEA.nCHN.lNsc)
summary(ESEA.nCHN.gomp)

ESEA.nCHN.RLGL2.AIC.vec <- c(AICc(ESEA.nCHN.rick), AICc(ESEA.nCHN.gomp))
ESEA.nCHN.RLGL2.dAIC.vec <- delta.IC(ESEA.nCHN.RLGL2.AIC.vec)
ESEA.nCHN.RLGL2.wAIC.vec <- weight.IC(ESEA.nCHN.RLGL2.dAIC.vec)
print(ESEA.nCHN.RLGL2.wAIC.vec)
ESEA.nCHN.ER.RLGL2 <- ESEA.nCHN.RLGL2.wAIC.vec[1]/ESEA.nCHN.RLGL2.wAIC.vec[2]
ESEA.nCHN.ER.RLGL2
1/ESEA.nCHN.ER.RLGL2

## calculate 95% prediction confidence intervals for Ricker and Gompertz
plot(ESEA.nCHN.use$Nsc, ESEA.nCHN.use$r, xlab="scaled N", ylab="r", pch=19)

# Ricker
ESEA.nCHN.x.pred.R <- data.frame(Nsc = seq(min(ESEA.nCHN.use$Nsc,na.rm=T), max(ESEA.nCHN.use$Nsc,na.rm=T), length.out=100))
ESEA.nCHN.pred.interval.R <- predict(ESEA.nCHN.rick, newdata=ESEA.nCHN.x.pred.R, interval="confidence", level=0.95)
lines(ESEA.nCHN.x.pred.R$Nsc, ESEA.nCHN.pred.interval.R[,"lwr"], col="blue", lty=2)
lines(ESEA.nCHN.x.pred.R$Nsc, ESEA.nCHN.pred.interval.R[,"upr"], col="blue", lty=2)

# Gompertz
ESEA.nCHN.x.pred.G <- data.frame(ESEA.nCHN.lNsc = seq(min(ESEA.nCHN.lNsc,na.rm=T), max(ESEA.nCHN.lNsc,na.rm=T), length.out=100))
ESEA.nCHN.pred.interval.G <- predict(ESEA.nCHN.gomp, newdata=ESEA.nCHN.x.pred.G, interval="confidence", level=0.95)
lines(exp(ESEA.nCHN.x.pred.G$ESEA.nCHN.lNsc), ESEA.nCHN.pred.interval.G[,"lwr"], col="green", lty=2)
lines(exp(ESEA.nCHN.x.pred.G$ESEA.nCHN.lNsc), ESEA.nCHN.pred.interval.G[,"upr"], col="green", lty=2)

# unscale Nsc
ESEA.nCHN.R.N.descaled <- attr(ESEA.nCHN.NSc, 'scaled:scale') * ESEA.nCHN.x.pred.R$Nsc
ESEA.nCHN.G.N.descaled <- attr(ESEA.nCHN.NSc, 'scaled:scale') * exp(ESEA.nCHN.x.pred.G$ESEA.nCHN.lNsc)

## export results
ESEA.nCHN.results.out <- data.frame(x.R=ESEA.nCHN.R.N.descaled,
                              R.fit=ESEA.nCHN.pred.interval.R[,"fit"],
                              R.upr=ESEA.nCHN.pred.interval.R[,"upr"],
                              R.lwr=ESEA.nCHN.pred.interval.R[,"lwr"],
                              x.G=ESEA.nCHN.G.N.descaled,
                              G.fit=ESEA.nCHN.pred.interval.G[,"fit"],
                              G.upr=ESEA.nCHN.pred.interval.G[,"upr"],
                              G.lwr=ESEA.nCHN.pred.interval.G[,"lwr"])
head(ESEA.nCHN.results.out)

# save to .csv
write.csv(ESEA.nCHN.results.out, "ESEA_minus_CHN_Ricker_Gompertz_fits.csv", row.names=F)




########################################################################
## linear models for examining contribution of per-capita consumption ##
########################################################################

## & population size to temperature anomaly
# import data
TAconN <- read.table("consump.csv", sep=",", header=T)
head(TAconN)
tail(TAconN)

# plot
par(mfrow=c(1,3))
plot(TAconN$pcEconsum, TAconN$TaMED, pch=19, xlab="per-capita consumption", ylab="temperature anomaly")
plot(TAconN$pop, TAconN$TaMED, pch=19, xlab="population size", ylab="temperature anomaly")
plot(TAconN$pop, TAconN$pcEconsum, pch=19, xlab="population size", ylab="per-capita consumption")
par(mfrow=c(1,1))

## Newey-West coefficient confidence intervals
pop_conN <- lm(TaMED ~ pop, data=TAconN)
summary(pop_conN)
pop_conN.coeftest <- coeftest(pop_conN, vcov = NeweyWest(pop_conN, prewhite = FALSE), save=T)
pop_conN.coeftest
pop_conN.coeftest.confint <- confint(pop_conN.coeftest)
pop_conN.coeftest.confint
pop_conN.coeftest.fitted <- attr(pop_conN.coeftest, "object")[[5]]
pop_conN.R2.NW <- cor(na.omit(TAconN$TaMED), pop_conN.coeftest.fitted)^2
pop_conN.R2.NW

## resampled
dat.use <- data.frame(pop=TAconN$pop,
                      anomLO=TAconN$TaLO,
                      anomUP=TAconN$TaUP)
head(dat.use)


iter <- 10000
itdiv <- iter/10

R2.vec <- ER.vec <- p.NW.vec <- slope.NW.vec <- R2.NW.vec <- rep(NA, iter)
for (i in 1:iter) {
  ta.samp <- runif(dim(dat.use)[1], dat.use$anomLO, dat.use$anomUP)
  R2.vec[i] <- linreg.ER(dat.use$pop, ta.samp)[2]
  ER.vec[i] <- linreg.ER(dat.use$pop, ta.samp)[1]
  
  popXta.fit <- lm(ta.samp ~ dat.use$pop)
  coeftest.fit <- coeftest(popXta.fit, vcov = NeweyWest(popXta.fit, prewhite = FALSE), save=T)
  
  coeftest.fitted <- attr(coeftest.fit, "object")[[5]]
  R2.NW.vec[i] <- cor(ta.samp, coeftest.fitted)^2
  slope.NW.vec[i] <- as.numeric(attr(coeftest.fit, "object")[[1]][2])
  p.NW.vec[i] <- coeftest.fit[7]
  
  if (i %% itdiv==0) print(i)
  
} # end i

R2lo <- quantile(R2.vec, probs=0.025, na.rm=T)
R2up <- quantile(R2.vec, probs=0.975, na.rm=T)
ERlo <- quantile(ER.vec, probs=0.025, na.rm=T)
ERup <- quantile(ER.vec, probs=0.975, na.rm=T)

print(c(R2lo, R2up))
print(c(ERlo, ERup))

p.NW.lo <- quantile(p.NW.vec, probs=0.025, na.rm=T)
p.NW.up <- quantile(p.NW.vec, probs=0.975, na.rm=T)
R2.NW.lo <- quantile(R2.NW.vec, probs=0.025, na.rm=T)
R2.NW.up <- quantile(R2.NW.vec, probs=0.975, na.rm=T)
slope.NW.lo <- quantile(slope.NW.vec, probs=0.025, na.rm=T)
slope.NW.up <- quantile(slope.NW.vec, probs=0.975, na.rm=T)

print(c(R2.NW.lo, R2.NW.up))
print(c(p.NW.lo, p.NW.up))
print(c(slope.NW.lo, slope.NW.up))



con_conN <- lm(TaMED ~ pcEconsum, data=TAconN)
summary(con_conN)
con_conN.coeftest <- coeftest(con_conN, vcov = NeweyWest(con_conN, prewhite = FALSE), save=T)
con_conN.coeftest
con_conN.coeftest.confint <- confint(con_conN.coeftest)
con_conN.coeftest.confint
con_conN.coeftest.fitted <- attr(con_conN.coeftest, "object")[[5]]
con_conN.R2.NW <- cor(na.omit(TAconN$TaMED), con_conN.coeftest.fitted)^2
con_conN.R2.NW


# models
mod1 <- "TaMED~pop+pcEconsum"
mod2 <- "TaMED~pop"
mod3 <- "TaMED~pcEconsum"
mod4 <- "TaMED~1"

## model vector
mod.vec <- c(mod1,mod2,mod3,mod4)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=TAconN, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit.sat <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=TAconN, na.action=na.omit)

check_model(fit.sat, detrend=F)
diagsat <- check_model(fit.sat, detrend=F)
diagsat$VIF$x
diagsat$VIF$y
plot_model(fit.sat, show.values=T, vline.color = "purple")


## resampling within temperature anomaly confidence interval
iter <- 10000
itdiv <- iter/10

mod1 <- "Ta.it~pop+con"
mod2 <- "Ta.it~pop"
mod3 <- "Ta.it~con"
mod4 <- "Ta.it~1"

## model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

topmod.vec <- wAICc.mod1.vec <- wAICc.mod2.vec <- wAICc.mod3.vec <- wAICc.mod4.vec <-
  DE.mod1.vec <- DE.mod2.vec <- DE.mod3.vec <- DE.mod4.vec <- rep(NA,iter)

for (j in 1:iter) {
  Ta.it <- runif(dim(TAconN)[1],min=TAconN$TaLO, max=TAconN$TaUP)
  TAconN.it1 <- data.frame(TAconN$pcEconsum, TAconN$pop, Ta.it)
  colnames(TAconN.it1) <- c("con", "pop", "Ta.it")
  
  # resample iterated dataset with replacement
  TAconN.it <- TAconN.it1[sample(1:30, replace=T), ]
  
  mod.list <- list()
  for(i in 1:length(mod.vec)) {
    fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=TAconN.it, na.action=na.omit)
    assign(paste("fit",i,sep=""), fit)
    mod.list[[i]] <- fit
  }
  sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
  row.names(sumtable) <- mod.vec
  summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
  summary.table
  
  topmod.vec[j] <- which(mod.vec == row.names(summary.table)[1])
                  
  wAICc.mod1.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[1]),5]
  wAICc.mod2.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[2]),5]
  wAICc.mod3.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[3]),5]
  wAICc.mod4.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[4]),5]
  
  DE.mod1.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[1]),9]
  DE.mod2.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[2]),9]
  DE.mod3.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[3]),9]
  DE.mod4.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[4]),9]
  
  if (j %% itdiv==0) print(j) 
  
} # end j

topmod.header <- mod.vec[as.numeric(attr(table(topmod.vec), "names"))]
topmod.table <- as.data.frame(table(topmod.vec)/iter)
topmod.table[,1] <- topmod.header
colnames(topmod.table) <- c("model", "%top")
topmod.table

mod1.AICmed <- median(wAICc.mod1.vec, na.rm=T)
mod2.AICmed <- median(wAICc.mod2.vec, na.rm=T)
mod3.AICmed <- median(wAICc.mod3.vec, na.rm=T)
mod4.AICmed <- median(wAICc.mod4.vec, na.rm=T)

mod1.AIClo <- quantile(wAICc.mod1.vec, probs=0.025, na.rm=T)
mod2.AIClo <- quantile(wAICc.mod2.vec, probs=0.025, na.rm=T)
mod3.AIClo <- quantile(wAICc.mod3.vec, probs=0.025, na.rm=T)
mod4.AIClo <- quantile(wAICc.mod4.vec, probs=0.025, na.rm=T)

mod1.AICup <- quantile(wAICc.mod1.vec, probs=0.975, na.rm=T)
mod2.AICup <- quantile(wAICc.mod2.vec, probs=0.975, na.rm=T)
mod3.AICup <- quantile(wAICc.mod3.vec, probs=0.975, na.rm=T)
mod4.AICup <- quantile(wAICc.mod4.vec, probs=0.975, na.rm=T)

mod1.AICmed.st <- mod1.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod2.AICmed.st <- mod2.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod3.AICmed.st <- mod3.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod4.AICmed.st <- mod4.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
AICmed.st.vec <- c(mod1.AICmed.st,mod2.AICmed.st,mod3.AICmed.st,mod4.AICmed.st)

mod1.AIClo.st <- mod1.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod2.AIClo.st <- mod2.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod3.AIClo.st <- mod3.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod4.AIClo.st <- mod4.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
AIClo.st.vec <- c(mod1.AIClo.st,mod2.AIClo.st,mod3.AIClo.st,mod4.AIClo.st)

mod1.AICup.st <- mod1.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod2.AICup.st <- mod2.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod3.AICup.st <- mod3.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod4.AICup.st <- mod4.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
AICup.st.vec <- c(mod1.AICup.st,mod2.AICup.st,mod3.AICup.st,mod4.AICup.st)

topmod.index <- c(which(topmod.table$model == mod.vec[1]),which(topmod.table$model == mod.vec[2]),
                  which(topmod.table$model == mod.vec[3]))
topmod.table$AICmdst <- round(AICmed.st.vec[topmod.index],4)
topmod.table$AIClost <- round(AIClo.st.vec[topmod.index],4)
topmod.table$AICupst <- round(AICup.st.vec[topmod.index],4)

mod1.DEmed <- median(DE.mod1.vec, na.rm=T)
mod2.DEmed <- median(DE.mod2.vec, na.rm=T)
mod3.DEmed <- median(DE.mod3.vec, na.rm=T)
mod4.DEmed <- median(DE.mod4.vec, na.rm=T)
DEmed.st.vec <- c(mod1.DEmed,mod2.DEmed,mod3.DEmed,mod4.DEmed)

topmod.table$DEmdst <- round(DEmed.st.vec[topmod.index], 1)
topmod.sort <- topmod.table[order(topmod.table[,2],decreasing=T),]
topmod.sort



## boosted regression tree (median TA)
brt.fit <- gbm.step(TAconN, gbm.x = attr(TAconN, "names")[c(2:3)], gbm.y = attr(TAconN, "names")[4],
                    family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001,
                    bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

brt.CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
brt.CV.cor.se <- 100 * brt.fit$cv.statistics$correlation.se
print(c(brt.CV.cor, brt.CV.cor.se))


# resampled BRT loop
biter <- 1000
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, 2, biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                             attr(TAconN, "names")[c(2:3)], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- N.ri <- E.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among years
  resamp.sub <- sort(sample(x = 1:dim(TAconN)[1], size = dim(TAconN)[1], replace=TRUE))
  dat.resamp <- TAconN[resamp.sub,]
  dat.resamp$TA.resamp <- runif(dim(dat.resamp)[1],min=dat.resamp$TaLO, max=dat.resamp$TaUP)
  
  # boosted regression tree
  brt.fit <- gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[c(2:3)], gbm.y = attr(dat.resamp, "names")[11],
                      family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75,
                      tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  
  length(summ.fit[[1]])
  
  if (length(summ.fit[[1]]) == 2) {
    # variable relative importance
    E.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][1])]
    N.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][2])]

    D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) /
                brt.fit$cv.statistics$deviance.mean
    D2.vec[b] <- D2
    CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
    CV.cor.vec[b] <- CV.cor
    CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
    CV.cor.se.vec[b] <- CV.cor.se
    
    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=2)
    ## output average predictions
    for (p in 1:2) {
      RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit$var.names
    
    val.arr[, , b] <- as.matrix(RESP.val.dat)
    pred.arr[, , b] <- as.matrix(RESP.pred.dat)
    
    print(b)
  }
  
  if (length(summ.fit[[1]]) != 2) {
    b <- b+1
    print(b)
  }
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | 
                                  pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(1,2)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) TA (higher→)", xlab="(←higher) E (lower→)")
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) TA (higher→)",  xlab="(←lower) N (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")
par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]
E.ri.update <- E.ri[1:biter]
N.ri.update <- N.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  E.mean <- mean(E.ri.update, na.rm=T); E.sd <- sd(E.ri.update, na.rm=T)
  N.mean <- mean(N.ri.update, na.rm=T); N.sd <- sd(N.ri.update, na.rm=T)

  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    E.ri.update[u] <- ifelse((E.ri.update[u] < (E.mean-kappa*E.sd) | E.ri.update[u] > (E.mean+kappa*E.sd)), NA, E.ri.update[u])
    N.ri.update[u] <- ifelse((N.ri.update[u] < (N.mean-kappa*N.sd) | N.ri.update[u] > (N.mean+kappa*N.sd)), NA, N.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

E.ri.lo <- quantile(E.ri.update, probs=0.025, na.rm=TRUE)
E.ri.med <- median(E.ri.update, na.rm=TRUE)
E.ri.up <- quantile(E.ri.update, probs=0.975, na.rm=TRUE)

N.ri.lo <- quantile(N.ri.update, probs=0.025, na.rm=TRUE)
N.ri.med <- median(N.ri.update, na.rm=TRUE)
N.ri.up <- quantile(N.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(E.ri.lo,N.ri.lo)
ri.med <- c(E.ri.med,N.ri.med)
ri.up <- c(E.ri.up,N.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(TAconN, "names")[c(2:3)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort


write.table(pred.med,file="NETA.BRT.boot.pred.med.csv",sep=",", row.names = T, col.names = T)
write.table(pred.lo,file="NETA.BRT.boot.pred.lo.csv",sep=",", row.names = T, col.names = T)
write.table(pred.up,file="NETA.BRT.boot.pred.up.csv",sep=",", row.names = T, col.names = T)
write.table(val.med,file="NETA.BRT.boot.val.med.csv",sep=",", row.names = T, col.names = T)
write.table(ri.sort,file="NETA.BRT.boot.ri.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.val.dat,file="NETA.BRT.val.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="NETA.BRT.pred.csv",sep=",", row.names = T, col.names = T)



#########################################
## & population size to global footprint
# import data
FPconN <- read.table("consump.csv", sep=",", header=T)
head(FPconN)
tail(FPconN)

# plot
par(mfrow=c(1,3))
plot(FPconN$pcEconsum, FPconN$footprint, pch=19, xlab="per-capita consumption", ylab="Earth's consumed")
plot(FPconN$pop, FPconN$footprint, pch=19, xlab="population size", ylab="Earth's consumed")
plot(FPconN$pop, FPconN$pcEconsum, pch=19, xlab="population size", ylab="per-capita consumption")
par(mfrow=c(1,1))

## Newey-West coefficient confidence intervals
pop_FPN <- lm(footprint ~ pop, data=FPconN)
summary(pop_FPN)
pop_FPN.coeftest <- coeftest(pop_FPN, vcov = NeweyWest(pop_FPN, prewhite = FALSE), save=T)
pop_FPN.coeftest
pop_FPN.coeftest.confint <- confint(pop_FPN.coeftest)
pop_FPN.coeftest.confint
pop_FPN.coeftest.fitted <- attr(pop_FPN.coeftest, "object")[[5]]
pop_FPN.R2.NW <- cor(na.omit(FPconN$footprint), pop_FPN.coeftest.fitted)^2
pop_FPN.R2.NW

## boosted regression tree (footprint)
brt.fit <- gbm.step(FPconN, gbm.x = attr(FPconN, "names")[c(2:3)], gbm.y = attr(FPconN, "names")[8],
                    family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001,
                    bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

brt.CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
brt.CV.cor.se <- 100 * brt.fit$cv.statistics$correlation.se
print(c(brt.CV.cor, brt.CV.cor.se))


# resampled BRT loop
biter <- 1000
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, 2, biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                           attr(TAconN, "names")[c(2:3)], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- N.ri <- E.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among years
  resamp.sub <- sort(sample(x = 1:dim(FPconN)[1], size = dim(FPconN)[1], replace=TRUE))
  dat.resamp <- FPconN[resamp.sub,]

  # boosted regression tree
  brt.fit <- gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[c(2:3)], gbm.y = attr(dat.resamp, "names")[8],
                      family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75,
                      tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  
  length(summ.fit[[1]])
  
  if (length(summ.fit[[1]]) == 2) {
    # variable relative importance
    E.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][1])]
    N.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][2])]
    
    D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) /
      brt.fit$cv.statistics$deviance.mean
    D2.vec[b] <- D2
    CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
    CV.cor.vec[b] <- CV.cor
    CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
    CV.cor.se.vec[b] <- CV.cor.se
    
    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=2)
    ## output average predictions
    for (p in 1:2) {
      RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit$var.names
    
    val.arr[, , b] <- as.matrix(RESP.val.dat)
    pred.arr[, , b] <- as.matrix(RESP.pred.dat)
    
    print(b)
  }
  
  if (length(summ.fit[[1]]) != 2) {
    b <- b+1
    print(b)
  }
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | 
                                  pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(1,2)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) FP (higher→)", xlab="(←higher) E (lower→)")
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) FP (higher→)",  xlab="(←lower) N (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")
par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]
E.ri.update <- E.ri[1:biter]
N.ri.update <- N.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  E.mean <- mean(E.ri.update, na.rm=T); E.sd <- sd(E.ri.update, na.rm=T)
  N.mean <- mean(N.ri.update, na.rm=T); N.sd <- sd(N.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    E.ri.update[u] <- ifelse((E.ri.update[u] < (E.mean-kappa*E.sd) | E.ri.update[u] > (E.mean+kappa*E.sd)), NA, E.ri.update[u])
    N.ri.update[u] <- ifelse((N.ri.update[u] < (N.mean-kappa*N.sd) | N.ri.update[u] > (N.mean+kappa*N.sd)), NA, N.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

E.ri.lo <- quantile(E.ri.update, probs=0.025, na.rm=TRUE)
E.ri.med <- median(E.ri.update, na.rm=TRUE)
E.ri.up <- quantile(E.ri.update, probs=0.975, na.rm=TRUE)

N.ri.lo <- quantile(N.ri.update, probs=0.025, na.rm=TRUE)
N.ri.med <- median(N.ri.update, na.rm=TRUE)
N.ri.up <- quantile(N.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(E.ri.lo,N.ri.lo)
ri.med <- c(E.ri.med,N.ri.med)
ri.up <- c(E.ri.up,N.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(FPconN, "names")[c(2:3)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort


write.table(pred.med,file="FP.BRT.boot.pred.med.csv",sep=",", row.names = T, col.names = T)
write.table(pred.lo,file="FP.BRT.boot.pred.lo.csv",sep=",", row.names = T, col.names = T)
write.table(pred.up,file="FP.BRT.boot.pred.up.csv",sep=",", row.names = T, col.names = T)
write.table(val.med,file="FP.BRT.boot.val.med.csv",sep=",", row.names = T, col.names = T)
write.table(ri.sort,file="FP.BRT.boot.ri.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.val.dat,file="FP.BRT.val.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="FP.BRT.pred.csv",sep=",", row.names = T, col.names = T)


########################################
## population and per-capita consumption relative to total CO2-e emissions (Gt)
## https://ourworldindata.org/co2-and-greenhouse-gas-emissions
EMconN <- read.table("consump.csv", sep=",", header=T)
head(EMconN)
tail(EMconN)

# plot
par(mfrow=c(1,3))
plot(EMconN$pcEconsum, EMconN$footprint, pch=19, xlab="per-capita consumption", ylab="total emissions")
plot(EMconN$pop, EMconN$footprint, pch=19, xlab="population size", ylab="total emissions")
plot(EMconN$pop, EMconN$pcEconsum, pch=19, xlab="population size", ylab="per-capita consumption")
par(mfrow=c(1,1))

## Newey-West coefficient confidence intervals
pop_EMN <- lm(CO2eGt ~ pop, data=EMconN)
summary(pop_EMN)
pop_EMN.coeftest <- coeftest(pop_EMN, vcov = NeweyWest(pop_EMN, prewhite = FALSE), save=T)
pop_EMN.coeftest
pop_EMN.coeftest.confint <- confint(pop_EMN.coeftest)
pop_EMN.coeftest.confint
pop_EMN.coeftest.fitted <- attr(pop_EMN.coeftest, "object")[[5]]
pop_EMN.R2.NW <- cor(na.omit(EMconN$footprint), pop_EMN.coeftest.fitted)^2
pop_EMN.R2.NW

## boosted regression tree (footprint)
brt.fit <- gbm.step(EMconN, gbm.x = attr(EMconN, "names")[c(2:3)], gbm.y = attr(EMconN, "names")[10],
                    family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001,
                    bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

brt.CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
brt.CV.cor.se <- 100 * brt.fit$cv.statistics$correlation.se
print(c(brt.CV.cor, brt.CV.cor.se))


# resampled BRT loop
biter <- 1000
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, 2, biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                                           attr(TAconN, "names")[c(2:3)], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- N.ri <- E.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among years
  resamp.sub <- sort(sample(x = 1:dim(EMconN)[1], size = dim(EMconN)[1], replace=TRUE))
  dat.resamp <- EMconN[resamp.sub,]
  
  # boosted regression tree
  brt.fit <- gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[c(2:3)], gbm.y = attr(dat.resamp, "names")[10],
                      family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75,
                      tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  
  length(summ.fit[[1]])
  
  if (length(summ.fit[[1]]) == 2) {
    # variable relative importance
    E.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][1])]
    N.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][2])]
    
    D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) /
      brt.fit$cv.statistics$deviance.mean
    D2.vec[b] <- D2
    CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
    CV.cor.vec[b] <- CV.cor
    CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
    CV.cor.se.vec[b] <- CV.cor.se
    
    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=2)
    ## output average predictions
    for (p in 1:2) {
      RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit$var.names
    
    val.arr[, , b] <- as.matrix(RESP.val.dat)
    pred.arr[, , b] <- as.matrix(RESP.pred.dat)
    
    print(b)
  }
  
  if (length(summ.fit[[1]]) != 2) {
    b <- b+1
    print(b)
  }
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | 
                                  pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(1,2)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) EM (higher→)", xlab="(←higher) E (lower→)")
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) EM (higher→)",  xlab="(←lower) N (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")
par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]
E.ri.update <- E.ri[1:biter]
N.ri.update <- N.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  E.mean <- mean(E.ri.update, na.rm=T); E.sd <- sd(E.ri.update, na.rm=T)
  N.mean <- mean(N.ri.update, na.rm=T); N.sd <- sd(N.ri.update, na.rm=T)
  
  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    E.ri.update[u] <- ifelse((E.ri.update[u] < (E.mean-kappa*E.sd) | E.ri.update[u] > (E.mean+kappa*E.sd)), NA, E.ri.update[u])
    N.ri.update[u] <- ifelse((N.ri.update[u] < (N.mean-kappa*N.sd) | N.ri.update[u] > (N.mean+kappa*N.sd)), NA, N.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

E.ri.lo <- quantile(E.ri.update, probs=0.025, na.rm=TRUE)
E.ri.med <- median(E.ri.update, na.rm=TRUE)
E.ri.up <- quantile(E.ri.update, probs=0.975, na.rm=TRUE)

N.ri.lo <- quantile(N.ri.update, probs=0.025, na.rm=TRUE)
N.ri.med <- median(N.ri.update, na.rm=TRUE)
N.ri.up <- quantile(N.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(E.ri.lo,N.ri.lo)
ri.med <- c(E.ri.med,N.ri.med)
ri.up <- c(E.ri.up,N.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(EMconN, "names")[c(2:3)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort


write.table(pred.med,file="EM.BRT.boot.pred.med.csv",sep=",", row.names = T, col.names = T)
write.table(pred.lo,file="EM.BRT.boot.pred.lo.csv",sep=",", row.names = T, col.names = T)
write.table(pred.up,file="EM.BRT.boot.pred.up.csv",sep=",", row.names = T, col.names = T)
write.table(val.med,file="EM.BRT.boot.val.med.csv",sep=",", row.names = T, col.names = T)
write.table(ri.sort,file="EM.BRT.boot.ri.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.val.dat,file="EM.BRT.val.csv",sep=",", row.names = T, col.names = T)
write.table(RESP.pred.dat,file="EM.BRT.pred.csv",sep=",", row.names = T, col.names = T)




# age structure # UN Population Division (World Population Prospects 2022); both sexes population by single age
# 1950-2021
astruct <- read.csv("agestructure1950_2021.csv", header=T)
head(astruct)

totpop <- 1000*apply(astruct[,-1], MARGIN=1, sum, na.rm=T)
age.vec <- seq(0.5,100.5,1)
age.mn <- rep(NA,dim(astruct)[1])
for (i in 1:dim(astruct)[1]) {
  age.mn[i] <- sum((1000*astruct[i, 2:102]) * age.vec)/totpop[i]
} # end a
plot(astruct$year, age.mn, type="l", xlab="year", ylab="mean age (years)")
plot(totpop, age.mn, type="l", ylab="mean age (years)", xlab="human population size")

pl15 <- (1000*apply(astruct[,2:16], MARGIN=1, sum, na.rm=T))/totpop

age.mn.out <- data.frame(astruct$year, totpop/10^9, age.mn, pl15)
colnames(age.mn.out) <- c("year", "Ntot", "ageMN", "young")

plot(age.mn.out$Ntot, age.mn.out$young, type="l", xlab="human population size", ylab="proportion < 15 years")
age.mn.out$year[which(age.mn.out$young == max(age.mn.out$young))]



########################################################################
## inequality relative to global environmental deterioration          ##
########################################################################

# Set WID_DATA_DIR to the directory containing the WID_data_XX.csv files.
# These series are WID's published regional aggregates, not averages of
# national Gini coefficients.

# download full dataset from 'https://wid.world/data/#:~:text=DOWNLOAD%20FULL-,DATASET,-NO%20INDICATOR%20SELECTED'
# set wid.dir to the directory containing the WID_data_XX.csv files (unzip downloaded file above)
wid.dir <- 'XXX'  # replace with your local path to WID_data_XX.csv files
#wid.dir <- Sys.getenv("WID_DATA_DIR")
if (!nzchar(wid.dir) || !dir.exists(wid.dir)) {
  stop("Set WID_DATA_DIR to the directory containing WID_data_XX.csv files.")
}

wid.scopes <- data.frame(
  scope=c("World", "Sub-Saharan Africa", "Latin America (WID)", "MENA (WID)"),
  code=c("WO", "XF", "XL", "XN"),
  alignment=c("Global series",
              "Matches the project's Sub-Saharan Africa region",
              "Published WID aggregate; not an exact LAC match",
              "Published WID aggregate; not an exact NAWA match"),
  stringsAsFactors=F
)

# Do not construct Oceania, ESEA, EUNA, or CSA from national Ginis: a
# population-weighted mean omits between-country inequality.
read.wid.gini <- function(wid.dir, scope, code, alignment) {
  file <- file.path(wid.dir, paste0("WID_data_", code, ".csv"))
  if (!file.exists(file)) {
    stop(paste("Missing WID input:", file))
  }

  dat <- read.csv(file, sep=";", stringsAsFactors=F)
  dat <- dat[dat$country == code &
               dat$variable %in% c("gdiincj992", "gptincj992") &
               dat$percentile == "p0p100" &
               as.character(dat$age) == "992" &
               dat$pop == "j", ]

  if (nrow(dat) == 0) {
    stop(paste("No adult equal-split WID Gini series found for", code))
  }
  if (any(duplicated(dat[,c("variable", "year")]))) {
    stop(paste("Duplicate WID Gini observations found for", code))
  }

  data.frame(scope=scope,
             wid.code=code,
             alignment=alignment,
             inequality.measure=ifelse(dat$variable == "gdiincj992",
                                        "Disposable-income Gini",
                                        "Pre-tax national-income Gini"),
             year=as.integer(dat$year),
             gini=as.numeric(dat$value),
             stringsAsFactors=F)
}

wid.series <- do.call(rbind, lapply(seq_len(nrow(wid.scopes)), function(i) {
  read.wid.gini(wid.dir=wid.dir,
                scope=wid.scopes$scope[i],
                code=wid.scopes$code[i],
                alignment=wid.scopes$alignment[i])
}))
wid.series <- wid.series[order(wid.series$scope,
                               wid.series$inequality.measure,
                               wid.series$year), ]

# Global environmental indicators complete from 1965-2023; WID
# Gini series provide common 1980-2023 analysis window.
data.dir.GH <- 'XXX' # replace with your local path to the directory containing consump.csv
wid.environment <- read.csv(file.path(data.dir.GH, "consump.csv"), header=T)
wid.environment <- data.frame(year=wid.environment$year,
                              pop=wid.environment$pop,
                              pcEconsum=wid.environment$pcEconsum,
                              TaMED=wid.environment$TaMED,
                              footprint=wid.environment$footprint,
                              CO2eGt=wid.environment$CO2eGt)
wid.dat <- wid.series[wid.series$year %in% wid.environment$year, ]
wid.environment.match <- wid.environment[match(wid.dat$year,
                                                wid.environment$year), -1]
wid.dat <- cbind(wid.dat, wid.environment.match)
head(wid.dat)
#wid.dat <- wid.dat[order(c(wid.dat$scope, wid.dat$inequality.measure, wid.dat$year)), ]



########################################################################
## regional energy use and ecological footprint                        ##
########################################################################

# Country energy use is sourced from the World Bank and ecological footprint
# and population from the Global Footprint Network. Regions follow WID's
# country classification so the predictors align with WID regional Gini series.
wid.regional.environment.file <- file.path(data.dir.GH,
                                           "WID_regional_energy_footprint.csv")
wid.regional.sources.file <- file.path(data.dir.GH,
                                       "WID_regional_energy_footprint_sources.csv")
wid.refresh.regional.environment <- identical(
  tolower(Sys.getenv("WID_REFRESH_REGIONAL_ENVIRONMENT", "false")), "true"
)

wid.countries <- read.csv(file.path(wid.dir, "WID_countries.csv"),
                          sep=";",
                          stringsAsFactors=F)
wid.region.membership <- rbind(
  data.frame(scope="Sub-Saharan Africa",
             isoa2=wid.countries$alpha2[
               wid.countries$region == "Africa" &
                 wid.countries$region2 != "North Africa"
             ]),
  data.frame(scope="Latin America (WID)",
             isoa2=wid.countries$alpha2[
               wid.countries$region == "Americas" &
                 wid.countries$region2 != "North America"
             ]),
  data.frame(scope="MENA (WID)",
             isoa2=wid.countries$alpha2[
               (wid.countries$region == "Africa" &
                  wid.countries$region2 == "North Africa") |
                 (wid.countries$region == "Asia" &
                    wid.countries$region2 == "West Asia")
             ])
)
wid.region.membership <- wid.region.membership[
  !is.na(wid.region.membership$isoa2) &
    nchar(wid.region.membership$isoa2) == 2, ]
wid.region.membership$country.name <- wid.countries$shortname[
  match(wid.region.membership$isoa2, wid.countries$alpha2)
]

wid.gfn.get <- function(record) {
  handle <- curl::new_handle()
  curl::handle_setheaders(
    handle,
    "Origin"="https://data.footprintnetwork.org",
    "Referer"="https://data.footprintnetwork.org/",
    "Accept"="application/json",
    "User-Agent"="R regional inequality analysis"
  )
  response <- curl::curl_fetch_memory(
    paste0("https://api.footprintnetwork.org/v1/data/all/all/", record),
    handle=handle
  )
  if (response$status_code != 200) {
    stop(paste("Global Footprint Network API request failed for", record,
               "with HTTP status", response$status_code))
  }
  jsonlite::fromJSON(rawToChar(response$content))
}

wid.owid.energy.get <- function() {
  handle <- curl::new_handle()
  curl::handle_setheaders(handle, "User-Agent"="R regional inequality analysis")
  response <- curl::curl_fetch_memory(
    "https://ourworldindata.org/grapher/per-capita-energy-use.csv",
    handle=handle
  )
  if (response$status_code != 200) {
    stop(paste("OWID per-capita energy-use download failed with HTTP status",
               response$status_code))
  }
  dat <- read.csv(text=rawToChar(response$content), stringsAsFactors=F)
  data.frame(country.name=dat$Entity,
             year=as.integer(dat$Year),
             energy.gj.per.capita=as.numeric(dat$Per.capita.energy.consumption) *
               0.0036,
             stringsAsFactors=F)
}

if (!file.exists(wid.regional.environment.file) ||
    wid.refresh.regional.environment) {
  if (!requireNamespace("curl", quietly=TRUE) ||
      !requireNamespace("jsonlite", quietly=TRUE)) {
    stop("Install the curl and jsonlite packages to download regional energy and footprint data.")
  }

  footprint <- wid.gfn.get("EFCtot")
  population <- wid.gfn.get("pop")
  energy <- wid.owid.energy.get()

  footprint <- footprint[footprint$isoa2 %in% wid.region.membership$isoa2 &
                           footprint$year >= 1980 &
                           footprint$year <= 2023, ]
  population <- population[population$isoa2 %in% wid.region.membership$isoa2 &
                             population$year >= 1980 &
                             population$year <= 2023, ]
  footprint.key <- paste(footprint$isoa2, footprint$year)
  population.index <- match(footprint.key,
                            paste(population$isoa2, population$year))
  energy.name.aliases <- c(
    "Côte d'Ivoire"="Cote d'Ivoire",
    "Tanzania, United Republic of"="Tanzania",
    "Congo, Democratic Republic of"="Democratic Republic of Congo",
    "Cabo Verde"="Cape Verde",
    "Venezuela, Bolivarian Republic of"="Venezuela",
    "Türkiye"="Turkey",
    "Libyan Arab Jamahiriya"="Libya",
    "Russian Federation"="Russia",
    "State of Palestine"="Palestine"
  )
  energy.country.name <- footprint$countryName
  alias.index <- match(energy.country.name, names(energy.name.aliases))
  energy.country.name[!is.na(alias.index)] <- energy.name.aliases[alias.index[!is.na(alias.index)]]
  energy.index <- match(paste(energy.country.name, footprint$year),
                        paste(energy$country.name, energy$year))
  regional.country.dat <- data.frame(
    isoa2=footprint$isoa2,
    year=footprint$year,
    footprint.total.gha=footprint$value,
    population=population$value[population.index],
    energy.gj.per.capita=energy$energy.gj.per.capita[energy.index],
    stringsAsFactors=F
  )
  regional.country.dat <- regional.country.dat[
    complete.cases(regional.country.dat), ]
  regional.country.dat$scope <- wid.region.membership$scope[
    match(regional.country.dat$isoa2, wid.region.membership$isoa2)
  ]

  wid.regional.environment <- do.call(rbind, lapply(
    split(regional.country.dat,
          list(regional.country.dat$scope, regional.country.dat$year),
          drop=T),
    function(x) {
      data.frame(
        scope=x$scope[1],
        year=x$year[1],
        n.countries=nrow(x),
        population=sum(x$population),
        energy.gj.per.capita=sum(x$energy.gj.per.capita * x$population) /
          sum(x$population),
        footprint.total.gha=sum(x$footprint.total.gha),
        footprint.gha.per.capita=sum(x$footprint.total.gha) / sum(x$population),
        stringsAsFactors=F
      )
    }
  ))
  wid.regional.environment <- wid.regional.environment[
    order(wid.regional.environment$scope, wid.regional.environment$year), ]
  write.csv(wid.regional.environment, wid.regional.environment.file, row.names=F)

  wid.regional.sources <- data.frame(
    dataset=c("Per-capita primary-energy consumption",
              "National Footprint and Biocapacity Accounts"),
    provider=c("Our World in Data",
               "Global Footprint Network"),
    url=c("https://ourworldindata.org/grapher/per-capita-energy-use.csv",
          "https://api.footprintnetwork.org/v1/data/all/all/EFCtot"),
    unit=c("GJ per person; converted from kWh",
           "Global hectares"),
    stringsAsFactors=F
  )
  write.csv(wid.regional.sources, wid.regional.sources.file, row.names=F)
} else {
  wid.regional.environment <- read.csv(wid.regional.environment.file,
                                       stringsAsFactors=F)
}


########################################################################
## regional greenhouse-gas emissions                                  ##
########################################################################

# Climate Watch reports territorial emissions in tonnes of CO2-equivalents
# by sector. Its annual country coverage starts in 1990, limiting these
# regional models to the 1990--2023 WID overlap.
wid.regional.emissions.file <- file.path(data.dir.GH,
                                         "WID_regional_emissions.csv")
wid.regional.emissions.sources.file <- file.path(
  data.dir.GH, "WID_regional_emissions_sources.csv"
)
wid.refresh.regional.emissions <- identical(
  tolower(Sys.getenv("WID_REFRESH_REGIONAL_EMISSIONS", "false")), "true"
)

wid.owid.ghg.get <- function() {
  handle <- curl::new_handle()
  curl::handle_setheaders(handle, "User-Agent"="R regional inequality analysis")
  response <- curl::curl_fetch_memory(
    "https://ourworldindata.org/grapher/greenhouse-gas-emissions-by-sector.csv",
    handle=handle
  )
  if (response$status_code != 200) {
    stop(paste("OWID greenhouse-gas download failed with HTTP status",
               response$status_code))
  }
  dat <- read.csv(text=rawToChar(response$content), stringsAsFactors=F)
  sector.columns <- setdiff(names(dat), c("Entity", "Code", "Year"))
  sector.values <- as.matrix(dat[sector.columns])
  emissions.co2e.tonnes <- rowSums(sector.values, na.rm=T)
  emissions.co2e.tonnes[rowSums(!is.na(sector.values)) == 0] <- NA
  data.frame(
    country.name=dat$Entity,
    year=as.integer(dat$Year),
    emissions.co2e.tonnes=emissions.co2e.tonnes,
    stringsAsFactors=F
  )
}

if (!file.exists(wid.regional.emissions.file) ||
    wid.refresh.regional.emissions) {
  if (!requireNamespace("curl", quietly=TRUE)) {
    stop("Install the curl package to download regional greenhouse-gas data.")
  }

  ghg.name.aliases <- c(
    "DR Congo"="Democratic Republic of Congo",
    "Cabo Verde"="Cape Verde",
    "Swaziland"="Eswatini",
    "Russian Federation"="Russia",
    "Syrian Arab Republic"="Syria",
    "Palestine"="Palestine"
  )
  ghg.country.name <- wid.region.membership$country.name
  alias.index <- match(ghg.country.name, names(ghg.name.aliases))
  ghg.country.name[!is.na(alias.index)] <- ghg.name.aliases[
    alias.index[!is.na(alias.index)]
  ]
  ghg <- wid.owid.ghg.get()
  ghg <- ghg[
    ghg$country.name %in% ghg.country.name &
      ghg$year >= 1990 & ghg$year <= 2023 &
      !is.na(ghg$emissions.co2e.tonnes),
  ]
  ghg$scope <- wid.region.membership$scope[
    match(ghg$country.name, ghg.country.name)
  ]
  wid.regional.emissions <- do.call(rbind, lapply(
    split(ghg, list(ghg$scope, ghg$year), drop=T),
    function(x) {
      data.frame(
        scope=x$scope[1],
        year=x$year[1],
        n.countries=nrow(x),
        emissions.co2e.tonnes=sum(x$emissions.co2e.tonnes),
        stringsAsFactors=F
      )
    }
  ))
  wid.regional.emissions$key <- paste(wid.regional.emissions$scope,
                                      wid.regional.emissions$year)
  environment.key <- paste(wid.regional.environment$scope,
                           wid.regional.environment$year)
  environment.index <- match(wid.regional.emissions$key, environment.key)
  wid.regional.emissions$population <- wid.regional.environment$population[
    environment.index
  ]
  wid.regional.emissions$emissions.co2e.tonnes.per.capita <-
    wid.regional.emissions$emissions.co2e.tonnes /
    wid.regional.emissions$population
  wid.regional.emissions$key <- NULL
  wid.regional.emissions <- wid.regional.emissions[
    order(wid.regional.emissions$scope, wid.regional.emissions$year), ]
  write.csv(wid.regional.emissions, wid.regional.emissions.file, row.names=F)

  wid.regional.emissions.sources <- data.frame(
    dataset="Greenhouse gas emissions by sector",
    provider="Climate Watch (via Our World in Data)",
    url="https://ourworldindata.org/grapher/greenhouse-gas-emissions-by-sector.csv",
    unit="Tonnes of CO2-equivalents; sectors summed",
    stringsAsFactors=F
  )
  write.csv(wid.regional.emissions.sources,
            wid.regional.emissions.sources.file,
            row.names=F)
} else {
  wid.regional.emissions <- read.csv(wid.regional.emissions.file,
                                     stringsAsFactors=F)
}

wid.coverage <- do.call(rbind, lapply(split(wid.dat,
                                            list(wid.dat$scope,
                                                 wid.dat$inequality.measure),
                                            drop=T),
                                    function(x) {
  data.frame(scope=x$scope[1],
             wid.code=x$wid.code[1],
             alignment=x$alignment[1],
             inequality.measure=x$inequality.measure[1],
             n.years=nrow(x),
             first.year=min(x$year),
             last.year=max(x$year),
             stringsAsFactors=F)
}))

wid.coefficient <- function(fit, method, scope, code, alignment,
                            inequality.measure, environmental.indicator) {
  nw <- coeftest(fit, vcov=NeweyWest(fit, lag=1, prewhite=F))
  data.frame(scope=scope,
             wid.code=code,
             alignment=alignment,
             inequality.measure=inequality.measure,
             environmental.indicator=environmental.indicator,
             method=method,
             n=nrow(model.frame(fit)),
             estimate=nw["environmental", "Estimate"],
             std.error=nw["environmental", "Std. Error"],
             statistic=nw["environmental", "t value"],
             p.value=nw["environmental", "Pr(>|t|)"],
             stringsAsFactors=F)
}

wid.association <- function(dat, environmental.indicator) {
  dat <- dat[order(dat$year), ]
  level.dat <- data.frame(gini=dat$gini,
                          environmental=dat[[environmental.indicator]],
                          year=dat$year)
  level.fit <- lm(gini ~ environmental + year, data=level.dat)

  diff.dat <- data.frame(year=dat$year[-1],
                         gini=diff(dat$gini),
                         environmental=diff(dat[[environmental.indicator]]),
                         consecutive=diff(dat$year) == 1)
  diff.dat <- diff.dat[diff.dat$consecutive, ]
  diff.fit <- lm(gini ~ environmental, data=diff.dat)

  rbind(wid.coefficient(level.fit,
                        "Level model controlling for linear time trend",
                        dat$scope[1],
                        dat$wid.code[1],
                        dat$alignment[1],
                        dat$inequality.measure[1],
                        environmental.indicator),
        wid.coefficient(diff.fit,
                        "First-difference model",
                        dat$scope[1],
                        dat$wid.code[1],
                        dat$alignment[1],
                        dat$inequality.measure[1],
                        environmental.indicator))
}

wid.associations <- do.call(rbind, lapply(split(wid.dat,
                                                list(wid.dat$scope,
                                                     wid.dat$inequality.measure),
                                                drop=T),
                                        function(x) {
  do.call(rbind, lapply(c("TaMED", "footprint", "CO2eGt"), function(y) {
    wid.association(x, y)
  }))
}))

write.csv(wid.series, "WID_inequality_series.csv", row.names=F)
write.csv(wid.coverage, "WID_inequality_coverage.csv", row.names=F)
write.csv(wid.associations, "WID_inequality_associations.csv", row.names=F)

wid.plot.variables <- data.frame(
  variable=c("TaMED", "footprint", "CO2eGt"),
  label=c("Global temperature anomaly",
          "Global ecological footprint",
          "Global CO2-e emissions (Gt)"),
  stringsAsFactors=F
)

wid.plot.dat <- do.call(rbind, lapply(seq_len(nrow(wid.plot.variables)), function(i) {
  data.frame(scope=wid.dat$scope,
             inequality.measure=wid.dat$inequality.measure,
             environmental.indicator=wid.plot.variables$label[i],
             environmental.value=wid.dat[[wid.plot.variables$variable[i]]],
             gini=wid.dat$gini,
             stringsAsFactors=F)
}))
wid.plot.dat$scope <- factor(wid.plot.dat$scope, levels=wid.scopes$scope)
wid.plot.dat$inequality.measure <- factor(wid.plot.dat$inequality.measure,
                                          levels=c("Disposable-income Gini",
                                                   "Pre-tax national-income Gini"))
wid.plot.dat$environmental.indicator <- factor(wid.plot.dat$environmental.indicator,
                                               levels=wid.plot.variables$label)

wid.bivariate.plot <- ggplot(wid.plot.dat,
                             aes(x=environmental.value, y=gini)) +
  geom_point(alpha=0.65, size=1.4) +
  geom_smooth(method="lm", formula=y ~ x, se=T, colour="firebrick") +
  facet_grid(inequality.measure + environmental.indicator ~ scope, scales="free") +
  labs(x=NULL,
       y="Gini coefficient",
       title="Bivariate relationships between inequality and global environmental indicators",
       subtitle="Published WID series, 1980-2023") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text.y=element_text(angle=0))
print(wid.bivariate.plot)
ggsave("WID_inequality_bivariate_plots.png",
       plot=wid.bivariate.plot,
       width=16,
       height=18,
       units="in",
       dpi=300)

wid.safe.filename <- function(x) {
  x <- tolower(gsub("[^[:alnum:]]+", "_", x))
  gsub("^_|_$", "", x)
}

wid.individual.plot.files <- character()
for (inequality.measure in levels(wid.plot.dat$inequality.measure)) {
  for (environmental.indicator in levels(wid.plot.dat$environmental.indicator)) {
    dat.plot <- wid.plot.dat[wid.plot.dat$inequality.measure == inequality.measure &
                             wid.plot.dat$environmental.indicator == environmental.indicator, ]

    wid.individual.plot <- ggplot(dat.plot,
                                  aes(x=environmental.value, y=gini)) +
      geom_point(alpha=0.65, size=1.8) +
      geom_smooth(method="lm", formula=y ~ x, se=T, colour="firebrick") +
      facet_wrap(~scope, scales="free_x") +
      labs(x=environmental.indicator,
           y="Gini coefficient",
           title=paste(inequality.measure, "and", environmental.indicator),
           subtitle="Published WID series, 1980-2023") +
      theme_bw() +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    print(wid.individual.plot)

    file <- paste0("WID_",
                   wid.safe.filename(inequality.measure),
                   "_by_",
                   wid.safe.filename(environmental.indicator),
                   ".png")
    ggsave(file,
           plot=wid.individual.plot,
           width=16,
           height=5,
           units="in",
           dpi=300)
    wid.individual.plot.files <- c(wid.individual.plot.files, file)
  }
}



########################################################################
## boosted regression trees including global income inequality         ##
########################################################################

# The default matches the existing BRT resampling procedure. Set
# WID_BRT_ITER to a smaller integer for exploratory runs.
wid.brt.biter <- as.integer(Sys.getenv("WID_BRT_ITER", "1000"))
wid.brt.points <- 100
if (is.na(wid.brt.biter) || wid.brt.biter < 2) {
  stop("WID_BRT_ITER must be an integer of at least 2.")
}

wid.kappa.vector <- function(x, kappa=2, kappa.n=5) {
  update <- x
  for (k in seq_len(kappa.n)) {
    centre <- mean(update, na.rm=T)
    spread <- sd(update, na.rm=T)
    update[update < centre-kappa*spread | update > centre+kappa*spread] <- NA
  }
  update
}

wid.kappa.curves <- function(x, kappa=2, kappa.n=5) {
  update <- x
  for (k in seq_len(kappa.n)) {
    centre <- apply(update, MARGIN=c(1,2), mean, na.rm=T)
    spread <- apply(update, MARGIN=c(1,2), sd, na.rm=T)
    for (b in seq_len(dim(update)[3])) {
      update[,,b] <- ifelse(update[,,b] < centre-kappa*spread |
                               update[,,b] > centre+kappa*spread,
                             NA,
                             update[,,b])
    }
  }
  update
}

wid.brt.summary <- function(x) {
  x <- wid.kappa.vector(x)
  c(lower=unname(quantile(x, probs=0.025, na.rm=T)),
    median=median(x, na.rm=T),
    upper=unname(quantile(x, probs=0.975, na.rm=T)))
}

run.wid.brt <- function(dat, response, inequality.measure,
                        scope="World",
                        predictors=c("pop", "pcEconsum", "gini"),
                        biter=wid.brt.biter, points=wid.brt.points) {
  dat <- dat[dat$scope == scope &
             dat$inequality.measure == inequality.measure, ]
  dat <- dat[,c(predictors, response)]
  colnames(dat)[4] <- "response"
  dat <- na.omit(dat)

  if (nrow(dat) < 10) {
    stop(paste("Insufficient complete WID BRT observations for", response,
               "and", inequality.measure))
  }

  initial.fit <- gbm.step(dat,
                          gbm.x=seq_along(predictors),
                          gbm.y=4,
                          family="gaussian",
                          max.trees=100000,
                          tolerance=0.0001,
                          learning.rate=0.0001,
                          bag.fraction=0.75,
                          tree.complexity=2,
                          silent=T,
                          tolerance.method="auto")

  curve.values <- curve.predictions <- array(
    data=NA,
    dim=c(points, length(predictors), biter),
    dimnames=list(NULL, predictors, NULL)
  )
  relative.influence <- matrix(NA,
                               nrow=biter,
                               ncol=length(predictors),
                               dimnames=list(NULL, predictors))
  d2 <- cv.correlation <- cv.correlation.se <- rep(NA, biter)

  for (b in seq_len(biter)) {
    dat.resamp <- dat[sample(seq_len(nrow(dat)), size=nrow(dat), replace=T), ]
    fit <- gbm.step(dat.resamp,
                    gbm.x=seq_along(predictors),
                    gbm.y=4,
                    family="gaussian",
                    max.trees=100000,
                    tolerance=0.0001,
                    learning.rate=0.001,
                    bag.fraction=0.75,
                    tree.complexity=2,
                    silent=T,
                    tolerance.method="auto")
    fit.summary <- summary(fit)
    relative.influence[b,] <- fit.summary$rel.inf[
      match(predictors, fit.summary$var)
    ]
    d2[b] <- 100 * (fit$cv.statistics$deviance.mean -
                    fit$self.statistics$mean.resid) /
      fit$cv.statistics$deviance.mean
    cv.correlation[b] <- 100 * fit$cv.statistics$correlation.mean
    cv.correlation.se[b] <- 100 * fit$cv.statistics$correlation.se

    for (p in seq_along(predictors)) {
      response.curve <- plot.gbm(fit,
                                 i.var=p,
                                 continuous.resolution=points,
                                 return.grid=T)
      curve.values[,p,b] <- response.curve[,1]
      curve.predictions[,p,b] <- response.curve[,2]
    }
    print(b)
  }

  curve.predictions <- wid.kappa.curves(curve.predictions)
  response.curves <- do.call(rbind, lapply(seq_along(predictors), function(p) {
    data.frame(predictor=predictors[p],
               value=apply(curve.values[,p,], 1, median, na.rm=T),
               fitted.lower=apply(curve.predictions[,p,], 1, quantile,
                                  probs=0.025, na.rm=T),
               fitted.median=apply(curve.predictions[,p,], 1, median, na.rm=T),
               fitted.upper=apply(curve.predictions[,p,], 1, quantile,
                                  probs=0.975, na.rm=T))
  }))

  relative.influence <- apply(relative.influence, 2, wid.brt.summary)
  relative.influence <- data.frame(predictor=colnames(relative.influence),
                                   lower=relative.influence["lower",],
                                   median=relative.influence["median",],
                                   upper=relative.influence["upper",],
                                   row.names=NULL)
  fit.metrics <- data.frame(
    metric=c("D2", "CV correlation", "CV correlation SE"),
    rbind(wid.brt.summary(d2),
          wid.brt.summary(cv.correlation),
          wid.brt.summary(cv.correlation.se)),
    row.names=NULL
  )

  list(initial.fit=initial.fit,
       response.curves=response.curves,
       relative.influence=relative.influence,
       fit.metrics=fit.metrics)
}

wid.brt.specifications <- expand.grid(
  response=c("TaMED", "footprint", "CO2eGt"),
  inequality.measure=c("Disposable-income Gini",
                       "Pre-tax national-income Gini"),
  stringsAsFactors=F
)
wid.brt.results <- vector("list", nrow(wid.brt.specifications))

for (i in seq_len(nrow(wid.brt.specifications))) {
  response <- wid.brt.specifications$response[i]
  inequality.measure <- wid.brt.specifications$inequality.measure[i]
  wid.brt.results[[i]] <- run.wid.brt(wid.dat,
                                      response=response,
                                      inequality.measure=inequality.measure)

  output.prefix <- paste0("WID_BRT_",
                          wid.safe.filename(inequality.measure),
                          "_",
                          wid.safe.filename(response))
  write.csv(wid.brt.results[[i]]$response.curves,
            paste0(output.prefix, "_response_curves.csv"),
            row.names=F)
  write.csv(wid.brt.results[[i]]$relative.influence,
            paste0(output.prefix, "_relative_influence.csv"),
            row.names=F)
  write.csv(wid.brt.results[[i]]$fit.metrics,
            paste0(output.prefix, "_fit_metrics.csv"),
            row.names=F)
}



########################################################################
## regional ecological-footprint boosted regression trees             ##
########################################################################

# Regional climate anomaly is not yet available. These models use regional
# total ecological footprint as the response.
wid.regional.gini <- wid.series[wid.series$scope %in%
                                  wid.regional.environment$scope, ]
wid.regional.gini.key <- paste(wid.regional.gini$scope,
                               wid.regional.gini$year)
wid.regional.environment.key <- paste(wid.regional.environment$scope,
                                      wid.regional.environment$year)
wid.regional.environment.index <- match(wid.regional.gini.key,
                                        wid.regional.environment.key)
wid.regional.brt.dat <- data.frame(
  scope=wid.regional.gini$scope,
  inequality.measure=wid.regional.gini$inequality.measure,
  year=wid.regional.gini$year,
  gini=wid.regional.gini$gini,
  population=wid.regional.environment$population[wid.regional.environment.index],
  energy.gj.per.capita=wid.regional.environment$energy.gj.per.capita[
    wid.regional.environment.index
  ],
  footprint.total.gha=wid.regional.environment$footprint.total.gha[
    wid.regional.environment.index
  ],
  stringsAsFactors=F
)
wid.regional.brt.dat <- wid.regional.brt.dat[
  complete.cases(wid.regional.brt.dat), ]

wid.regional.brt.specifications <- expand.grid(
  scope=sort(unique(wid.regional.brt.dat$scope)),
  inequality.measure=c("Disposable-income Gini",
                       "Pre-tax national-income Gini"),
  stringsAsFactors=F
)
wid.regional.brt.results <- vector("list", nrow(wid.regional.brt.specifications))

for (i in seq_len(nrow(wid.regional.brt.specifications))) {
  scope <- wid.regional.brt.specifications$scope[i]
  inequality.measure <- wid.regional.brt.specifications$inequality.measure[i]
  wid.regional.brt.results[[i]] <- run.wid.brt(
    wid.regional.brt.dat,
    response="footprint.total.gha",
    inequality.measure=inequality.measure,
    scope=scope,
    predictors=c("population", "energy.gj.per.capita", "gini")
  )

  output.prefix <- paste0("WID_regional_footprint_BRT_",
                          wid.safe.filename(scope),
                          "_",
                          wid.safe.filename(inequality.measure))
  write.csv(wid.regional.brt.results[[i]]$response.curves,
            paste0(output.prefix, "_response_curves.csv"),
            row.names=F)
  write.csv(wid.regional.brt.results[[i]]$relative.influence,
            paste0(output.prefix, "_relative_influence.csv"),
            row.names=F)
  write.csv(wid.regional.brt.results[[i]]$fit.metrics,
            paste0(output.prefix, "_fit_metrics.csv"),
            row.names=F)
}



########################################################################
## regional CO2-e boosted regression trees                             ##
########################################################################

wid.regional.emissions.gini <- wid.series[
  wid.series$scope %in% wid.regional.emissions$scope, ]
wid.regional.emissions.gini.key <- paste(wid.regional.emissions.gini$scope,
                                         wid.regional.emissions.gini$year)
wid.regional.emissions.key <- paste(wid.regional.emissions$scope,
                                    wid.regional.emissions$year)
wid.regional.emissions.index <- match(wid.regional.emissions.gini.key,
                                      wid.regional.emissions.key)
wid.regional.emissions.environment.key <- paste(
  wid.regional.environment$scope, wid.regional.environment$year
)
wid.regional.emissions.environment.index <- match(
  wid.regional.emissions.gini.key, wid.regional.emissions.environment.key
)
wid.regional.emissions.brt.dat <- data.frame(
  scope=wid.regional.emissions.gini$scope,
  inequality.measure=wid.regional.emissions.gini$inequality.measure,
  year=wid.regional.emissions.gini$year,
  gini=wid.regional.emissions.gini$gini,
  population=wid.regional.environment$population[
    wid.regional.emissions.environment.index
  ],
  energy.gj.per.capita=wid.regional.environment$energy.gj.per.capita[
    wid.regional.emissions.environment.index
  ],
  emissions.co2e.tonnes=wid.regional.emissions$emissions.co2e.tonnes[
    wid.regional.emissions.index
  ],
  stringsAsFactors=F
)
wid.regional.emissions.brt.dat <- wid.regional.emissions.brt.dat[
  complete.cases(wid.regional.emissions.brt.dat), ]

wid.regional.emissions.brt.specifications <- expand.grid(
  scope=sort(unique(wid.regional.emissions.brt.dat$scope)),
  inequality.measure=c("Disposable-income Gini",
                       "Pre-tax national-income Gini"),
  stringsAsFactors=F
)
wid.regional.emissions.brt.results <- vector(
  "list", nrow(wid.regional.emissions.brt.specifications)
)

for (i in seq_len(nrow(wid.regional.emissions.brt.specifications))) {
  scope <- wid.regional.emissions.brt.specifications$scope[i]
  inequality.measure <- wid.regional.emissions.brt.specifications$inequality.measure[i]
  wid.regional.emissions.brt.results[[i]] <- run.wid.brt(
    wid.regional.emissions.brt.dat,
    response="emissions.co2e.tonnes",
    inequality.measure=inequality.measure,
    scope=scope,
    predictors=c("population", "energy.gj.per.capita", "gini")
  )

  output.prefix <- paste0("WID_regional_emissions_BRT_",
                          wid.safe.filename(scope),
                          "_",
                          wid.safe.filename(inequality.measure))
  write.csv(wid.regional.emissions.brt.results[[i]]$response.curves,
            paste0(output.prefix, "_response_curves.csv"),
            row.names=F)
  write.csv(wid.regional.emissions.brt.results[[i]]$relative.influence,
            paste0(output.prefix, "_relative_influence.csv"),
            row.names=F)
  write.csv(wid.regional.emissions.brt.results[[i]]$fit.metrics,
            paste0(output.prefix, "_fit_metrics.csv"),
            row.names=F)
}
