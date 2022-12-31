# Clear data and load packages
rm(list = ls())
library(survival)
library(boot)

set.seed(10122)

# Load data
dat <- read.csv("/Volumes/Kawut_CCDIPH/Surrogate/surrogate.csv",
                header = TRUE)

# Interpolate 12 and 24 week data to get missing 16 week data for surrogates
# Not done for COMPERA 2 or FPHR (non-inv) because Trial 1010 has no data
dat$reveal2_wk16[is.na(dat$reveal2_wk16)] <-
  (2*dat$reveal2_wk12[is.na(dat$reveal2_wk16)] +
     dat$reveal2_wk24[is.na(dat$reveal2_wk16)]) / 3
dat$reveal_lite_wk16[is.na(dat$reveal_lite_wk16)] <-
  (2*dat$reveal_lite_wk12[is.na(dat$reveal_lite_wk16)] +
     dat$reveal_lite_wk24[is.na(dat$reveal_lite_wk16)]) / 3
dat$compera_full_wk16[is.na(dat$compera_full_wk16)] <-
  (2*dat$compera_full_wk12[is.na(dat$compera_full_wk16)] +
     dat$compera_full_wk24[is.na(dat$compera_full_wk16)]) / 3
dat$compera_2_wk16[is.na(dat$compera_2_wk16)] <-
  (2*dat$compera_2_wk12[is.na(dat$compera_2_wk16)] +
     dat$compera_2_wk24[is.na(dat$compera_2_wk16)]) / 3

# Make numeric treatment variable
dat$trt <- 1*(dat$study_arm_assignment == "Active")

# Clean up clinical outcome variables
dat$cw_day_full <- dat$cw_day
dat$cw_day_full[is.na(dat$cw_day_full)] <- dat$esp_day[is.na(dat$cw_day_full)]
dat$cw_day_full[dat$cw_day_full == 0] <- 1

dat$death_day_full <- dat$death_day
dat$death_day_full[is.na(dat$death_day_full)] <-
  dat$esp_day[is.na(dat$death_day_full)]
dat$death_day_full[dat$death_day_full == 0] <- 1

dat$cw_bin <- 0
dat$cw_bin[dat$cw == "CLINICAL WORSENING #1"] <- 1
dat$death_bin <- 0
dat$death_bin[dat$death == "DEATH"] <- 1

# Make control variables new types
dat$studyid_f <- as.factor(dat$studyid)
dat$eti_code_f <- as.factor(dat$eti_code)
dat$sysbp <- as.numeric(dat$sysbp)

# Make scores new types
dat$compera_cat0 <- factor(round(dat$compera_full_wk0))
levels(dat$compera_cat0) <- c("Low", "Intermediate", "High")
dat$compera_cat16 <- factor(round(dat$compera_full_wk16))
levels(dat$compera_cat16) <- c("Low", "Intermediate", "High")

dat$compera2_cat0 <- factor(round(dat$compera_2_wk0))
levels(dat$compera2_cat0) <- c("Low", "Intermediate-Low",
                               "Intermediate-High", "High")
dat$compera2_cat16 <- factor(round(dat$compera_2_wk16))
levels(dat$compera2_cat16) <- c("Low", "Intermediate-Low",
                                "Intermediate-High", "High")

dat$fphr_cat0 <- 3 - dat$fphr_noninv_wk0
dat$fphr_cat0[dat$fphr_cat0 == 0] <- 1
dat$fphr_cat0 <- factor(dat$fphr_cat0)
levels(dat$fphr_cat0) <- c("Low", "Intermediate", "High")

dat$fphr_cat16 <- 3 - dat$fphr_noninv_wk16
dat$fphr_cat16[dat$fphr_cat16 == 0] <- 1
dat$fphr_cat16 <- factor(dat$fphr_cat16)
levels(dat$fphr_cat16) <- c("Low", "Intermediate", "High")



######################################################
# Functions
#####################################################

# point estimates
mediate_surv = function(model.mar,model.con,treat,data) {
  m.mar=coxph(model.mar,data=data)
  m.con=coxph(model.con,data=data)
  a=m.mar$coefficients[treat]
  b=m.con$coefficients[treat]
  NIE = a - b
  NDE = b
  TE  = a
  MP  = NIE/TE
  c(NIE=NIE,NDE=NDE,TE=TE,MP=MP)
}

# bootstrap confidence interval
# Input: same to mediate_surv. Here R is the number of bootstrap replicates (
#        which I suggest >=500, ideally 1000 or more).
mediate_surv_ci = function(model.mar,model.con,treat,data,R=1000) {
  ConstructBootFun = function(d,i) {
    mediate_surv(model.mar=model.mar,model.con=model.con,treat=treat,data=d[i,])
  }
  boot_res=suppressWarnings(boot(data, ConstructBootFun, R = R, stype = "i"))
  out <- as.data.frame(matrix(NA,ncol=3,nrow=4))
  out[,1]=boot_res$t0
  rownames(out)=names(boot_res$t0)
  for (j in (1:4)) {
    out[j,2]  = sd(boot_res$t[,j])
    out[j,3:4]=boot.ci(boot_res,type="perc",index=j)$percent[4:5]
  }
  colnames(out) <- c("Estimate","S.E","CI.lower","CI.upper")
  out
}

###########################################################################
# REVEAL 2.0 Mediation analyses
dat$reveal2_low <- 1*(round(dat$reveal2_wk16) <= 6)

# primary outcome
model.mar = Surv(cw_day_full, cw_bin) ~ trt + reveal2_wk0 + factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + reveal2_low + reveal2_wk0 + factor(studyid_f)
treat = "trt"

mediate_surv(model.mar,model.con,treat,dat)
mediate_surv_ci(model.mar,model.con,treat,dat,R=1000) # ideally R=1000

# secondary outcome
model.mar = Surv(death_day_full, death_bin) ~ trt + reveal2_wk0+factor(studyid_f)
model.con = Surv(death_day_full, death_bin) ~ trt + reveal2_low +reveal2_wk0 + factor(studyid_f)

#model.mar = Surv(death_day_full, death_bin) ~ trt + reveal2_wk0
#model.con = Surv(death_day_full, death_bin) ~ trt + reveal2_wk16+reveal2_wk0

treat = "trt"

mediate_surv(model.mar,model.con,treat,dat)
mediate_surv_ci(model.mar,model.con,treat,dat,R=100) # ideally R=1000



###########################################################################
# REVEAL Lite Mediation analyses
dat$reveal_lite_low <- 1*(round(dat$reveal_lite_wk16) <= 5)

# primary outcome
model.mar = Surv(cw_day_full, cw_bin) ~ trt + reveal_lite_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + reveal_lite_low +reveal_lite_wk0 + factor(studyid_f)

mediate_surv(model.mar,model.con,treat,dat)
mediate_surv_ci(model.mar,model.con,treat,dat,R=1000)



###########################################################################
# COMPERA Mediation analyses
dat$compera_full_low <- 1*(dat$compera_cat16 == "Low")

# primary outcome
model.mar = Surv(cw_day_full, cw_bin) ~ trt + compera_full_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + compera_full_low +compera_full_wk0 + factor(studyid_f)

mediate_surv(model.mar,model.con,treat,dat)
mediate_surv_ci(model.mar,model.con,treat,dat,R=1000)


###########################################################################
# COMPERA 2 Mediation analyses
dat$compera_2_low <- 1*(dat$compera2_cat16 == "Low")

# primary outcome
model.mar = Surv(cw_day_full, cw_bin) ~ trt + compera_2_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + compera_2_low +compera_2_wk0 + factor(studyid_f)

mediate_surv(model.mar,model.con,treat,dat)
mediate_surv_ci(model.mar,model.con,treat,dat,R=1000)



###########################################################################
# FPHR Mediation analyses
dat$fphr_noninv_low <- 1*(dat$fphr_cat16 == "Low")

# primary outcome
model.mar = Surv(cw_day_full, cw_bin) ~ trt + fphr_noninv_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + fphr_noninv_low +fphr_noninv_wk0 + factor(studyid_f)

mediate_surv(model.mar,model.con,treat,dat)
mediate_surv_ci(model.mar,model.con,treat,dat,R=1000)


######
# Treatment on survival alone
mod <- coxph(Surv(death_day_full, death_bin) ~ trt, data = dat)
exp(confint(mod))





############################################################################
# Sensitivity analysis: Use AFT models

# Functions
# point estimates
mediate_aft = function(model.mar,model.con,treat,data) {
  m.mar=survreg(model.mar,data=data,dist="weibull")
  m.con=survreg(model.con,data=data,dist="weibull")
  a=m.mar$coefficients[treat]
  b=m.con$coefficients[treat]
  NIE = a - b
  NDE = b
  TE  = a
  MP  = NIE/TE
  c(NIE=NIE,NDE=NDE,TE=TE,MP=MP)
}

# bootstrap confidence interval
# Input: same to mediate_aft. Here R is the number of bootstrap replicates (
#        which I suggest >=500, ideally 1000 or more).
mediate_aft_ci = function(model.mar,model.con,treat,data,R=1000) {
  ConstructBootFun = function(d,i) {
    mediate_aft(model.mar=model.mar,model.con=model.con,treat=treat,data=d[i,])
  }
  boot_res=suppressWarnings(boot(data, ConstructBootFun, R = R, stype = "i"))
  out <- as.data.frame(matrix(NA,ncol=3,nrow=4))
  out[,1]=boot_res$t0
  rownames(out)=names(boot_res$t0)
  for (j in (1:4)) {
    out[j,2]  = sd(boot_res$t[,j])
    out[j,3:4]=boot.ci(boot_res,type="perc",index=j)$percent[4:5]
  }
  colnames(out) <- c("Estimate","S.E","CI.lower","CI.upper")
  out
}

###########################################################################
# REVEAL 2.0 Mediation analyses
model.mar = Surv(cw_day_full, cw_bin) ~ trt + reveal2_wk0 + factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + reveal2_low + reveal2_wk0 + factor(studyid_f)

mediate_aft(model.mar,model.con,treat,dat)
mediate_aft_ci(model.mar,model.con,treat,dat,R=1000)


###########################################################################
# REVEAL Lite Mediation analyses
model.mar = Surv(cw_day_full, cw_bin) ~ trt + reveal_lite_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + reveal_lite_low +reveal_lite_wk0 + factor(studyid_f)

mediate_aft(model.mar,model.con,treat,dat)
mediate_aft_ci(model.mar,model.con,treat,dat,R=1000)



###########################################################################
# COMPERA Mediation analyses
model.mar = Surv(cw_day_full, cw_bin) ~ trt + compera_full_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + compera_full_low +compera_full_wk0 + factor(studyid_f)

mediate_aft(model.mar,model.con,treat,dat)
mediate_aft_ci(model.mar,model.con,treat,dat,R=1000)


###########################################################################
# COMPERA 2 Mediation analyses
model.mar = Surv(cw_day_full, cw_bin) ~ trt + compera_2_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + compera_2_low +compera_2_wk0 + factor(studyid_f)

mediate_aft(model.mar,model.con,treat,dat)
mediate_aft_ci(model.mar,model.con,treat,dat,R=1000)


###########################################################################
# FPHR Mediation analyses
model.mar = Surv(cw_day_full, cw_bin) ~ trt + fphr_noninv_wk0+factor(studyid_f)
model.con = Surv(cw_day_full, cw_bin) ~ trt + fphr_noninv_low +fphr_noninv_wk0 + factor(studyid_f)

mediate_aft(model.mar,model.con,treat,dat)
mediate_aft_ci(model.mar,model.con,treat,dat,R=1000)
