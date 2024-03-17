rm(list=ls())

# Required sources
library(INLA)
library(INLAjoint)
library(KMsurv)      # Dataset for AFT and PH examples
library(smcure)      # Dataset for Mixture Cure example
library(compeir)     # Dataset for Competing Risks example
library(p3state.msm) # Dataset for Multi-state example
library(frailtyHL)   # Dataset for Frailty example
library(JMbayes)     # Dataset for Joint Model 1 example
library(frailtypack) # Dataset for Joint Model 2 example


###############################################################################################
#                             1 - ACCELERATED FAILURE TIME MODELS                             #
###############################################################################################
data(larynx)
larynx$age <- as.numeric(scale(larynx$age))
larynx$diagyr <- as.numeric(scale(larynx$diagyr))
larynx$stage <- as.factor(larynx$stage)

# Fit the model
m1.weib <- joint(formSurv = inla.surv(time = time, event = delta) ~ stage + age + diagyr,
                 basRisk = "weibullsurv", dataSurv = larynx, control = list(config=TRUE))

# Posterior summary
m1.gumb <- sapply(m1.weib$marginals.fixed, function(m)
                  inla.zmarginal(inla.tmarginal(function(x) -x, m)))
gumb.scale <- inla.tmarginal(function(x) 1/x, m1.weib$marginals.hyperpar[[1]])
rbind("Gumbel scale"=inla.zmarginal(gumb.scale), t(m1.gumb))

# Posterior marginal distribution plots
par(mfrow=c(2,4), mar = c(3,3,3,1))
m1.marg <- append(list("Gumbel scale" = gumb.scale),
                  lapply(m1.weib$marginals.fixed, function(m) inla.tmarginal(function(x) -x, m)))
sapply(1:7, function (x) plot(m1.marg[[x]], type = "l", xlab = "", ylab = "",
                              main = paste0("Density of ", names(m1.marg)[x])))

# Posterior sample
m1.sample <- inla.posterior.sample(1, m1.weib)

# Non-linear transformation of marginal
inla.zmarginal(inla.tmarginal(function(x) exp(-x), m1.weib$marginals.fixed[["stage4_S1"]]))

###############################################################################################
#                               2 - PROPORTIONAL HAZARDS MODELS                               #
###############################################################################################
# Fit the model
m2.pwc <- joint(formSurv = inla.surv(time = time, event = delta) ~ stage + age + diagyr,
                basRisk = "rw1", NbasRisk = 3, dataSurv = larynx)

# Posterior summary
summary(m2.pwc)
summary(m2.pwc)$BaselineValues
summary(m2.pwc, hr=TRUE)

###############################################################################################
#                        2bis - STRATIFIED PROPORTIONAL HAZARDS MODELS                        #
###############################################################################################
larynx$stage <- as.integer(larynx$stage)
m2b.pwc <- joint(formSurv = inla.surv(time = time, event = delta) ~ age + diagyr,
                 dataSurv = larynx, control=list(strata=list("stage")))

# Posterior summary
plot(m2b.pwc)
# makes sense to have higher baseline risk for each stage of cancer as it indicates how bad the disease is.

###############################################################################################
#                                   3 - MIXTURE CURE MODELS                                   #
###############################################################################################
data(bmt)
# Fit the model
m3.cure <- joint(formSurv = inla.surv(time = Time, event = Status,
                                      cure = cbind("Int"=1, "TRT"=TRT)) ~ TRT,
                 basRisk = "weibullsurv", dataSurv = bmt)

# Posterior summary
summary(m3.cure)

# compute cure fraction for the two treatment arms
smp.cure <- inla.hyperpar.sample(500, m3.cure)[, 2:3]
quantile(inv.logit(smp.cure[,1]), c(0.025, 0.5, 0.975)) # cure fraction allogenic
quantile(inv.logit(rowSums(smp.cure)), c(0.025, 0.5, 0.975)) # cure fraction autologous

# HR of treatment for non-cured
summary(m3.cure, hr=T)$SurvEff[[1]]["TRT_S1",]

###############################################################################################
#                                 4 - COMPETING RISKS MODELS                                  #
###############################################################################################
data(okiss)
delta <- matrix(c(as.integer(okiss$status == 1), as.integer(okiss$status == 2),
                  as.integer(okiss$status == 7)), ncol = 3)
head(delta)

# Fit the model
m4.cr <- joint(formSurv = list(inla.surv(time = time, event = delta[,1]) ~ allo + sex,
                               inla.surv(time = time, event = delta[,2]) ~ allo + sex,
                               inla.surv(time = time, event = delta[,3]) ~ allo + sex),
               basRisk = rep("weibullsurv", 3), dataSurv = okiss)

# Posterior summary
summary(m4.cr)
summary(m4.cr, hr=T)$SurvEff[[1]]["allo_S1",]

# Cumulative Incidence Functions
riskW <- function(t, lambda, alpha) lambda*alpha*t^(alpha-1)
t <- seq(0, 100, len=100)
risk1 <- riskW(t, exp(m4.cr$summary.fixed["Intercept_S1", "mean"]), m4.cr$summary.hyperpar$mean[1])
risk2 <- riskW(t, exp(m4.cr$summary.fixed["Intercept_S2", "mean"]), m4.cr$summary.hyperpar$mean[2])
risk3 <- riskW(t, exp(m4.cr$summary.fixed["Intercept_S3", "mean"]), m4.cr$summary.hyperpar$mean[3])
surv <- exp(-cumsum(risk1)-cumsum(risk2)-cumsum(risk3))
CIF <- cbind(cumsum(risk1*surv), cumsum(risk2*surv), cumsum(risk3*surv))

plot(t, CIF[,1], type="l", ylim=c(0,1), lwd=3, main="Cumulative Incidence Functions", ylab="Probability of event", xlab="Time")
lines(t, CIF[,2], lty=2, col=2, lwd=3)
lines(t, CIF[,3], lty=3, col=3, lwd=3)
legend("top", c("infection", "death", "end of neutropenia"),y.intersp = 0.5, x.intersp = 0.5, lty=c(1,3,2), lwd=c(3,3,3), col=c(1,3,2), horiz=TRUE)

###############################################################################################
#                                   5 - MULTI-STATE MODELS                                    #
###############################################################################################
data(heart2)
event <- matrix(c(heart2$delta, heart2$status * (1 - heart2$delta),
                  heart2$delta * heart2$status), ncol = 3)
head(event)

# Fit the model
heart2$times3 <- ifelse(heart2$times2 == 0, heart2$times1, heart2$times2)
m5.ms <- joint(formSurv = list(inla.surv(times1, event[,1]) ~ age + year + surgery,
                               inla.surv(times3, event[,2]) ~ age + year + surgery,
                               inla.surv(times2, event[,3]) ~ age + year + surgery),
               basRisk = rep("weibullsurv", 3), dataSurv = heart2)

# Posterior summary
summary(m5.ms)

# No clock-reset specification
# surv.obj.ms3 <- inla.surv(time = heart2$time, truncation = heart2$times1, event = event[,3])

# probability transitions
t <- seq(0.1, 1000, by=1)
riskW <- function(t, lambda, alpha) lambda*alpha*t^(alpha-1)
risk1 <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S1", "mean"]),
               m5.ms$summary.hyperpar$mean[1])
risk2 <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S2", "mean"]),
               m5.ms$summary.hyperpar$mean[2])
risk3 <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S3", "mean"]),
               m5.ms$summary.hyperpar$mean[3])
risk1s <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S1", "mean"] +
                       m5.ms$summary.fixed["surgery_S1", "mean"]),
                m5.ms$summary.hyperpar$mean[1])
risk2s <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S2", "mean"] +
                       m5.ms$summary.fixed["surgery_S2", "mean"]),
                m5.ms$summary.hyperpar$mean[2])
risk3s <- riskW(t, exp(m5.ms$summary.fixed["Intercept_S3", "mean"] +
                       m5.ms$summary.fixed["surgery_S3", "mean"]),
                m5.ms$summary.hyperpar$mean[3])
p11 <- exp(-cumsum(risk1)-cumsum(risk2))
p11s <- exp(-cumsum(risk1s)-cumsum(risk2s))
p22 <- exp(-cumsum(risk3))
p22s <- exp(-cumsum(risk3s))
p12 <- cumsum(p11*risk1*p22)
p12s <- cumsum(p11s*risk1s*p22s)
p13 <- 1-p11-p12
p13s <- 1-p11s-p12s
p23 <- 1-p22
p23s <- 1-p22s

par(mfrow=c(2,3), oma = c(1.5, 1.5, 1, 1), mai = c(0.2, 0.2, 0.3, 0.1), xpd=NA)
plot(t, p11, type="l", ylim=c(0,1), lwd=2, main="Probability of 1->1", xlab="", ylab="")
lines(t, p11s, col=2, lty=2, lwd=2)
plot(t, p12, type="l", ylim=c(0,1), lwd=2, main="Probability of 1->2", xlab="", ylab="")
lines(t, p12s, col=2, lty=2, lwd=2)
plot(t, p13, type="l", ylim=c(0,1), lwd=2, main="Probability of 1->3", xlab="", ylab="")
lines(t, p13s, col=2, lty=2, lwd=2)
plot(t, p22, type="l", ylim=c(0,1), lwd=2, main="Probability of 2->2", xlab="", ylab="")
lines(t, p22s, col=2, lty=2, lwd=2)
plot(t, p23, type="l", ylim=c(0,1), lwd=2, main="Probability of 2->3", xlab="", ylab="")
lines(t, p23s, col=2, lty=2, lwd=2)
legend(1400, 0.6, c("No prior surgery", "Prior surgery"), lty=c(1,2), lwd=c(2,2), col=c(1,2))

###############################################################################################
#                                     6 - FRAILTY MODELS                                      #
###############################################################################################
library(frailtyHL)
data(kidney)
kidney$sex <- kidney$sex - 1   # 0: male (reference)

# Fit the model
m6.frlt <- joint(formSurv = inla.surv(time = time, event = status) ~ sex + (1 | id),
                 basRisk = "weibullsurv", id = "id", dataSurv = kidney)

# Posterior summary
summary(m6.frlt)
summary(m6.frlt, hr=T)$SurvEff[[1]]["sex_S1",]

###############################################################################################
#                      7 - JOINT MODELS OF LONGITUDINAL AND SURVIVAL DATA                     #
###############################################################################################
data(prothro)
data(prothros)

# Fit the model
m7.jm1 <- joint(formSurv = inla.surv(time = Time, event = death) ~ treat, dataSurv = prothros,
                formLong = pro ~ time + treat + (1 + time | id), dataLong = prothro,
                basRisk = "weibullsurv", family = "lognormal", id = "id",
                timeVar = "time", assoc = "SRE")

# Posterior summary
summary(m7.jm1, sdcor = TRUE)   # Standard deviation instead of variance

# Posterior marginal distribution plots
m7.plots <- plot(m7.jm1, sdcor = FALSE)   # Variance instead of standard deviation
# Longitudinal fixed effects and variance error term
m7.plots$Outcomes$L1
# Survival fixed effect
m7.plots$Outcomes$S1
# Association parameter
m7.plots$Associations
# Variance-covariance matrix of the random effects
m7.plots$Covariances
# Weibull baseline hazard parameters
m7.plots$Baseline

# Multiplicative effect of treatment on prothrombin levels
inla.zmarginal(inla.tmarginal(function(x) exp(x), m7.jm1$marginals.fixed[["treatprednisone_L1"]]))

###############################################################################################
#  8 - JOINT MODELS OF LONGITUDINAL SEMICONTINUOUS, RECURRENT EVENTS AND TERMINAL EVENT DATA  #
###############################################################################################
data(colorectal)
data(colorectalLongi)
colorectalLongi$y <- round((colorectalLongi$tumor.size*0.3+1)^(1/0.3), 5)
colorectalLongiPositive <- colorectalLongi[colorectalLongi$y > 0,]
colorectalLongi$z <- ifelse(colorectalLongi$y == 0, 0, 1)
colorectalSurv <- subset(colorectal, new.lesions == 0)

# Fit the model
m8.jm2 <- joint(formSurv = list(inla.surv(time = time1, truncation = time0,  event = new.lesions) ~ treatment + (1 | id),
                                inla.surv(time = time1, event = state) ~ treatment),
                                formLong = list(z ~ year * treatment + (1 | id),
                                                y ~ year * treatment + (1 + year | id)),
                                dataSurv = list(colorectal, colorectalSurv),
                                dataLong = list(colorectalLongi, colorectalLongiPositive),
                                basRisk = c("rw2","rw2"), assocSurv = TRUE, timeVar = "year",
                                family = c("binomial","lognormal"), id = "id", corLong = TRUE,
                                assoc = list(c("SRE_ind","SRE_ind"),c("SRE_ind","SRE_ind")),
                control=list(int.strategy="eb"))

# Posterior summary
summary(m8.jm2)

# Posterior baseline hazard curves in log10 scale
plot(m8.jm2)$Baseline + scale_y_log10()

# Clock-reset specification
# surv.obj.jm2.recu <- inla.surv(time = colorectal$gap.time, event = colorectal$new.lesions)

# Probability of positive tumor size 3 years after randomization for an average patient with sequential treatment
inv.logit(m8.jm2$summary.fixed["Intercept_L1", 1]+3*m8.jm2$summary.fixed["year_L1", 1])

# Probability of positive tumor size 3 years after randomization for an average patient with combination treatment
inv.logit(m8.jm2$summary.fixed["Intercept_L1", 1]+3*m8.jm2$summary.fixed["year_L1", 1]+
          m8.jm2$summary.fixed["treatmentC_L1", 1]+3*m8.jm2$summary.fixed["year.X.treatmentC_L1", 1])

# Multiplicative effect of treatment on the average tumor size evolution over time, conditional on observing a positive tumor size
inla.zmarginal(inla.tmarginal(function(x) exp(x), m8.jm2$marginals.fixed[[10]]))














