rm(list=ls())

# Required sources
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
surv.obj <- inla.surv(time = larynx$time, event = larynx$delta)
m1.weib <- joint(formSurv = surv.obj ~ stage + age + diagyr, basRisk = "weibullsurv",
                            dataSurv = larynx, control = list(variant=0, cfg=T))

# Posterior summary
m1.gumb <- sapply(m1.weib$marginals.fixed, function(m) inla.zmarginal(inla.tmarginal(function(x) -x, m), silent = T))
post.m1.gumb <- t(m1.gumb[-c(4,6),])
colnames(post.m1.gumb) <- colnames(m1.weib$summary.hyperpar[-6])
rbind(m1.weib$summary.hyperpar[-6], post.m1.gumb)

# Posterior marginal distribution plots
par(mfrow=c(2,4))
m1.marg <- append(lapply(m1.weib$marginals.fixed, function(m) inla.tmarginal(function(x) -x, m)),
                   list("Gumbel scale" = inla.tmarginal(function(x) 1/x, m1.weib$marginals.hyperpar[[1]])))
sapply(1:7, function (x) plot(m1.marg[[x]], type = "l", xlab = "", ylab = "", main = paste0("Density of ", names(m1.marg)[x])))

# Posterior sample
m1.sample <- inla.posterior.sample(1, m1.weib)


###############################################################################################
#                               2 - PROPORTIONAL HAZARDS MODELS                               #
###############################################################################################
# Fit the model
surv.obj <- inla.surv(time = larynx$time, event = larynx$delta)
m2.pwc <- joint(formSurv = surv.obj ~ stage + age + diagyr, basRisk = "rw1",
                NbasRisk = 3, dataSurv = larynx)

# Posterior summary
summary(m2.pwc)
summary(m2.pwc)$BaselineValues


###############################################################################################
#                                   3 - MIXTURE CURE MODELS                                   #
###############################################################################################
data(bmt)

# Fit the model
surv.obj.cure <- inla.surv(time = bmt$Time, event = bmt$Status,
                           cure = cbind("Int"=1, "TRT"=bmt$TRT))
m3.cure <- joint(formSurv = surv.obj.cure ~ TRT, basRisk = "weibullsurv",
                 dataSurv = bmt, control = list(variant=0))

# Posterior summary
summary(m3.cure)


###############################################################################################
#                                 4 - COMPETING RISKS MODELS                                  #
###############################################################################################
data(okiss)
delta <- matrix(c(as.integer(okiss$status == 1), as.integer(okiss$status == 2),
                  as.integer(okiss$status == 7)), ncol = 3)
head(delta)

# Fit the model
surv.obj.cr1 <- inla.surv(time = okiss$time, event = delta[,1])
surv.obj.cr2 <- inla.surv(time = okiss$time, event = delta[,2])
surv.obj.cr3 <- inla.surv(time = okiss$time, event = delta[,3])
m4.cr <- joint(formSurv = list(surv.obj.cr1 ~ allo + sex,
                               surv.obj.cr2 ~ allo + sex,
                               surv.obj.cr3 ~ allo + sex),
               basRisk = c("weibullsurv", "weibullsurv", "weibullsurv"),
               dataSurv = okiss, control = list(variant=0))

# Posterior summary
summary(m4.cr)


###############################################################################################
#                                   5 - MULTI-STATE MODELS                                    #
###############################################################################################
data(heart2)
event <- matrix(c(heart2$delta, heart2$status * (1 - heart2$delta),
                  heart2$delta * heart2$status), ncol = 3)
head(event)

# Fit the model
heart2$times3 <- ifelse(heart2$times2 == 0, heart2$times1, heart2$times2)
surv.obj.ms1 <- inla.surv(time = heart2$times1, event = event[,1])
surv.obj.ms2 <- inla.surv(time = heart2$times3, event = event[,2])
surv.obj.ms3 <- inla.surv(time = heart2$times2, event = event[,3])
m5.ms <- joint(formSurv = list(surv.obj.ms1 ~ age + year + surgery,
                               surv.obj.ms2 ~ age + year + surgery,
                               surv.obj.ms3 ~ age + year + surgery),
               basRisk = c("weibullsurv", "weibullsurv", "weibullsurv"),
               dataSurv = heart2, control = list(variant=0))

# Posterior summary
summary(m5.ms)

# No clock-reset specification
# surv.obj.ms3 <- inla.surv(time = heart2$time, truncation = heart2$times1, event = event[,3])


###############################################################################################
#                                     6 - FRAILTY MODELS                                      #
###############################################################################################
library(frailtyHL)
data(kidney)
kidney$sex <- kidney$sex - 1   # 0: male (reference)

# Fit the model
surv.obj.frlt <- inla.surv(time = kidney$time, event = kidney$status)
m6.frlt <- joint(formSurv = surv.obj.frlt ~ sex + (1 | id), basRisk = "weibullsurv",
                 id = "id", dataSurv = kidney, control = list(variant=0))

# Posterior summary
summary(m6.frlt)


###############################################################################################
#                      7 - JOINT MODELS OF LONGITUDINAL AND SURVIVAL DATA                     #
###############################################################################################
data(prothro)
data(prothros)

# Fit the model
surv.obj.jm1 <- inla.surv(time = prothros$Time, event = prothros$death)
m7.jm1 <- joint(formSurv = surv.obj.jm1 ~ treat, dataSurv = prothros,
                formLong = pro ~ time + treat + (1 + time | id), dataLong = prothro,
                basRisk = "weibullsurv", family = "lognormal", id = "id",
                timeVar = "time", assoc = "SRE", control = list(variant=0))

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
surv.obj.jm2.recu <- inla.surv(time = colorectal$time1, truncation = colorectal$time0,
                               event = colorectal$new.lesions)
surv.obj.jm2.term <- inla.surv(time = colorectalSurv$time1, event = colorectalSurv$state)
m8.jm2 <- joint(formSurv = list(surv.obj.jm2.recu ~ treatment + (1 | id),
                                surv.obj.jm2.term ~ treatment),
                                formLong = list(z ~ year * treatment + (1 | id),
                                                y ~ year * treatment + (1 + year | id)),
                                dataSurv = list(colorectal, colorectalSurv),
                                dataLong = list(colorectalLongi, colorectalLongiPositive),
                                basRisk = c("rw2","rw2"), assocSurv = TRUE, timeVar = "year",
                                family = c("binomial","lognormal"), id = "id", corLong = TRUE,
                                assoc = list(c("SRE_ind","SRE_ind"),c("SRE_ind","SRE_ind")))

# Posterior summary
summary(m8.jm2)

# Posterior baseline hazard curves in log10 scale
plot(m8.jm2)$Baseline + scale_y_log10()

# Clock-reset specification
# surv.obj.jm2.recu <- inla.surv(time = colorectal$gap.time, event = colorectal$new.lesions)
