# Clear data and load libraries
rm(list = ls())
library(dplyr)
library(flextable)
library(ggplot2)
library(officer)
library(survival)
library(MASS)

# Load data
dat <- read.csv("/Volumes/Kawut_CCDIPH/Surrogate/surrogate.csv",
                header = TRUE)

# Interpolate 12 and 24 week data to get missing 16 week data for surrogates
# Not done for FPHR (non-inv) because Trial 1010 has no data
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

dat$cw_bin <- 0
dat$cw_bin[dat$cw == "CLINICAL WORSENING #1"] <- 1
dat$death_bin <- 0
dat$death_bin[dat$death == "DEATH"] <- 1

# Create region variable
r1list <- c("Canada", "United States")
r2list <- c("Australia", "Austria", "Belgium", "France", "Germany", "Greece",
            "Italy", "Japan", "Netherlands", "Spain", "Sweden",
            "United Kingdom")
r3list <- c("Argentina", "Canada", "Chile", "Columbia", "Mexico", "Peru",
            "U.S.A")
r4list <- c("Australia", "Austria", "Belarus", "Belgium", "Bulgaria",
            "Croatia", "Denmark", "France", "Germany", "Hungary", "Israel",
            "Italy", "Netherlands", "Norway", "Poland", "Romania", "Russia",
            "Serbia", "Slovakia", "South Africa", "Spain", "Sweden", "Turkey",
            "U.K.", "Ukraine")
r5list <- c("China", "HongKong", "India", "Malaysia", "Singapore", "Taiwan",
            "Thailand")
r6list <- c("ARG", "CAN", "CHL", "COL", "MEX", "PER", "USA")
r7list <- c("AUS", "AUT", "BEL", "BLR", "CHE", "CZE", "DEU", "DNK", "ESP",
            "FRA", "GBR", "GRC", "HUN", "IRL", "ISR", "ITA", "NLD", "POL",
            "ROU", "RUS", "SRB", "SVK", "SWE", "TUR", "UKR")
r8list <- c("CHN", "IND", "KOR", "MYS", "SGP", "THA", "TWN")

dat$studyregion <- ""
dat$studyregion[dat$studyid == 1004 & dat$country %in% r1list] <-
  "Trial 1: North America"
dat$studyregion[dat$studyid == 1004 & dat$country %in% r2list] <-
  "Trial 1: Europe/Australia"
dat$studyregion[dat$studyid == 1010 & dat$country %in% r3list] <-
  "Trial 2: Americas"
dat$studyregion[dat$studyid == 1010 & dat$country %in% r4list] <-
  "Trial 2: Europe/Australia"
dat$studyregion[dat$studyid == 1010 & dat$country %in% r5list] <-
  "Trial 2: Asia"
dat$studyregion[dat$studyid == 1013 & dat$country %in% r6list] <-
  "Trial 3: Americas"
dat$studyregion[dat$studyid == 1013 & dat$country %in% r7list] <-
  "Trial 3: Europe/Australia"
dat$studyregion[dat$studyid == 1013 & dat$country %in% r8list] <-
  "Trial 3: Asia"
dat$studyregion <- as.factor(dat$studyregion)

# Subset to analysis data set
dat <- dat[dat$cw_day_full >= 7*16, ]


###########################################################################
# REVEAL 2.0 Meta-analyses
# Trial-region specific effects and associations
dat$reveal2_low <- 1*(round(dat$reveal2_wk16) <= 6)

coef_rev_med <- list()
coef_rev_out <- list()
for (i in 1:8) {

  med_mod <- glm(reveal2_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])

  coef_rev_med[[i]] <- summary(med_mod)$coefficients
  coef_rev_out[[i]] <- summary(out_mod)$coefficients

}

trialdat_rev <-
  data.frame(surreff = c(coef_rev_med[[1]][2, 1], coef_rev_med[[2]][2, 1],
                         coef_rev_med[[3]][2, 1], coef_rev_med[[4]][2, 1],
                         coef_rev_med[[5]][2, 1], coef_rev_med[[6]][2, 1],
                         coef_rev_med[[7]][2, 1], coef_rev_med[[8]][2, 1]),
             clineff = c(coef_rev_out[[1]][1, 1], coef_rev_out[[2]][1, 1],
                         coef_rev_out[[3]][1, 1], coef_rev_out[[4]][1, 1],
                         coef_rev_out[[5]][1, 1], coef_rev_out[[6]][1, 1],
                         coef_rev_out[[7]][1, 1], coef_rev_out[[8]][1, 1]),
             w = c(1/coef_rev_out[[1]][1, 3]^2, 1/coef_rev_out[[2]][1, 3]^2,
                   1/coef_rev_out[[3]][1, 3]^2, 1/coef_rev_out[[4]][1, 3]^2,
                   1/coef_rev_out[[5]][1, 3]^2, 1/coef_rev_out[[6]][1, 3]^2,
                   1/coef_rev_out[[7]][1, 3]^2, 1/coef_rev_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_rev <- lm(clineff ~ surreff, data = trialdat_rev, weights = w)
r2label_rev <- round(summary(trialmod_rev)$r.squared, 2)

# Leave-one-out analysis
loo_preds_rev <- rep(NA, 8)
for (i in 1:8) {

  trialmod_rev <- lm(clineff ~ surreff, data = trialdat_rev[-i, ], weights = w)
  loo_preds_rev[i] <- predict(trialmod_rev, newdata = trialdat_rev[i, ])

}


###########################################################################
# REVEAL Lite Meta-analyses
# Trial-region specific effects and associations
dat$reveal_lite_low <- 1*(round(dat$reveal_lite_wk16) <= 5)

coef_revlite_med <- list()
coef_revlite_out <- list()
for (i in 1:8) {

  med_mod <- glm(reveal_lite_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])

  coef_revlite_med[[i]] <- summary(med_mod)$coefficients
  coef_revlite_out[[i]] <- summary(out_mod)$coefficients

}

trialdat_revlite <-
  data.frame(surreff = c(coef_revlite_med[[1]][2, 1],
                         coef_revlite_med[[2]][2, 1],
                         coef_revlite_med[[3]][2, 1],
                         coef_revlite_med[[4]][2, 1],
                         coef_revlite_med[[5]][2, 1],
                         coef_revlite_med[[6]][2, 1],
                         coef_revlite_med[[7]][2, 1],
                         coef_revlite_med[[8]][2, 1]),
             clineff = c(coef_revlite_out[[1]][1, 1],
                         coef_revlite_out[[2]][1, 1],
                         coef_revlite_out[[3]][1, 1],
                         coef_revlite_out[[4]][1, 1],
                         coef_revlite_out[[5]][1, 1],
                         coef_revlite_out[[6]][1, 1],
                         coef_revlite_out[[7]][1, 1],
                         coef_revlite_out[[8]][1, 1]),
             w = c(1/coef_revlite_out[[1]][1, 3]^2,
                   1/coef_revlite_out[[2]][1, 3]^2,
                   1/coef_revlite_out[[3]][1, 3]^2,
                   1/coef_revlite_out[[4]][1, 3]^2,
                   1/coef_revlite_out[[5]][1, 3]^2,
                   1/coef_revlite_out[[6]][1, 3]^2,
                   1/coef_revlite_out[[7]][1, 3]^2,
                   1/coef_revlite_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_revlite <-
  lm(clineff ~ surreff, data = trialdat_revlite, weights = w)
r2label_revlite <- round(summary(trialmod_revlite)$r.squared, 2)

# Leave-one-out analysis
loo_preds_revlite <- rep(NA, 8)
for (i in 1:8) {

  trialmod_revlite <- lm(clineff ~ surreff,
                         data = trialdat_revlite[-i, ], weights = w)
  loo_preds_revlite[i] <- predict(trialmod_revlite,
                                  newdata = trialdat_revlite[i, ])

}


###########################################################################
# COMPERA Meta-analyses
# Trial-region specific effects and associations
dat$compera_cat16 <- factor(round(dat$compera_full_wk16 + 0.01))
levels(dat$compera_cat16) <- c("Low", "Intermediate", "High")
dat$compera_full_low <- 1*(dat$compera_cat16 == "Low")

coef_comp_med <- list()
coef_comp_out <- list()
for (i in 1:8) {

  med_mod <- glm(compera_full_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])

  coef_comp_med[[i]] <- summary(med_mod)$coefficients
  coef_comp_out[[i]] <- summary(out_mod)$coefficients

}

trialdat_comp <-
  data.frame(surreff = c(coef_comp_med[[1]][2, 1], coef_comp_med[[2]][2, 1],
                         coef_comp_med[[3]][2, 1], coef_comp_med[[4]][2, 1],
                         coef_comp_med[[5]][2, 1], coef_comp_med[[6]][2, 1],
                         coef_comp_med[[7]][2, 1], coef_comp_med[[8]][2, 1]),
             clineff = c(coef_comp_out[[1]][1, 1], coef_comp_out[[2]][1, 1],
                         coef_comp_out[[3]][1, 1], coef_comp_out[[4]][1, 1],
                         coef_comp_out[[5]][1, 1], coef_comp_out[[6]][1, 1],
                         coef_comp_out[[7]][1, 1], coef_comp_out[[8]][1, 1]),
             w = c(1/coef_comp_out[[1]][1, 3]^2, 1/coef_comp_out[[2]][1, 3]^2,
                   1/coef_comp_out[[3]][1, 3]^2, 1/coef_comp_out[[4]][1, 3]^2,
                   1/coef_comp_out[[5]][1, 3]^2, 1/coef_comp_out[[6]][1, 3]^2,
                   1/coef_comp_out[[7]][1, 3]^2, 1/coef_comp_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_comp <- lm(clineff ~ surreff, data = trialdat_comp, weights = w)
r2label_comp <- round(summary(trialmod_comp)$r.squared, 2)

# Leave-one-out analysis
loo_preds_comp <- rep(NA, 8)
for (i in 1:8) {

  trialmod_comp <- lm(clineff ~ surreff,
                      data = trialdat_comp[-i, ], weights = w)
  loo_preds_comp[i] <- predict(trialmod_comp,
                                  newdata = trialdat_comp[i, ])

}


###########################################################################
# COMPERA 2.0 Meta-analyses
# Trial-region specific effects and associations
dat$compera2_cat16 <- factor(round(dat$compera_2_wk16 + 0.01))
levels(dat$compera2_cat16) <- c("Low", "Intermediate-Low",
                                "Intermediate-High", "High")
dat$compera_2_low <- 1*(dat$compera2_cat16 == "Low")

coef_comp2_med <- list()
coef_comp2_out <- list()
for (i in 1:8) {

  med_mod <- glm(compera_2_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])

  coef_comp2_med[[i]] <- summary(med_mod)$coefficients
  coef_comp2_out[[i]] <- summary(out_mod)$coefficients

}

trialdat_comp2 <-
  data.frame(surreff = c(coef_comp2_med[[1]][2, 1], coef_comp2_med[[2]][2, 1],
                         coef_comp2_med[[3]][2, 1], coef_comp2_med[[4]][2, 1],
                         coef_comp2_med[[5]][2, 1], coef_comp2_med[[6]][2, 1],
                         coef_comp2_med[[7]][2, 1], coef_comp2_med[[8]][2, 1]),
             clineff = c(coef_comp2_out[[1]][1, 1], coef_comp2_out[[2]][1, 1],
                         coef_comp2_out[[3]][1, 1], coef_comp2_out[[4]][1, 1],
                         coef_comp2_out[[5]][1, 1], coef_comp2_out[[6]][1, 1],
                         coef_comp2_out[[7]][1, 1], coef_comp2_out[[8]][1, 1]),
             w = c(1/coef_comp2_out[[1]][1, 3]^2,
                   1/coef_comp2_out[[2]][1, 3]^2,
                   1/coef_comp2_out[[3]][1, 3]^2,
                   1/coef_comp2_out[[4]][1, 3]^2,
                   1/coef_comp2_out[[5]][1, 3]^2,
                   1/coef_comp2_out[[6]][1, 3]^2,
                   1/coef_comp2_out[[7]][1, 3]^2,
                   1/coef_comp2_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_comp2 <- lm(clineff ~ surreff, data = trialdat_comp2, weights = w)
r2label_comp2 <- round(summary(trialmod_comp2)$r.squared, 2)

# Leave-one-out analysis
loo_preds_comp2 <- rep(NA, 8)
for (i in 1:8) {

  trialmod_comp2 <- lm(clineff ~ surreff,
                      data = trialdat_comp2[-i, ], weights = w)
  loo_preds_comp2[i] <- predict(trialmod_comp2,
                               newdata = trialdat_comp2[i, ])

}


###########################################################################
# FPHR Meta-analyses
# Trial-region specific effects and associations
dat$fphr_cat16 <- 3 - dat$fphr_noninv_wk16
dat$fphr_cat16[dat$fphr_cat16 == 0] <- 1
dat$fphr_cat16 <- factor(dat$fphr_cat16)
levels(dat$fphr_cat16) <- c("Low", "Intermediate", "High")
dat$fphr_noninv_low <- 1*(dat$fphr_cat16 == "Low")

coef_fphr_med <- list()
coef_fphr_out <- list()
for (i in c(1, 2, 6, 7, 8)) {

  med_mod <- glm(fphr_noninv_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])

  coef_fphr_med[[i]] <- summary(med_mod)$coefficients
  coef_fphr_out[[i]] <- summary(out_mod)$coefficients

}

trialdat_fphr <-
  data.frame(surreff = c(coef_fphr_med[[1]][2, 1], coef_fphr_med[[2]][2, 1],
                         coef_fphr_med[[6]][2, 1], coef_fphr_med[[7]][2, 1],
                         coef_fphr_med[[8]][2, 1]),
             clineff = c(coef_fphr_out[[1]][1, 1], coef_fphr_out[[2]][1, 1],
                         coef_fphr_out[[6]][1, 1], coef_fphr_out[[7]][1, 1],
                         coef_fphr_out[[8]][1, 1]),
             w = c(1/coef_fphr_out[[1]][1, 3]^2,
                   1/coef_fphr_out[[2]][1, 3]^2,
                   1/coef_fphr_out[[6]][1, 3]^2,
                   1/coef_fphr_out[[7]][1, 3]^2,
                   1/coef_fphr_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)[c(1, 2, 6, 7, 8)]))

trialmod_fphr <- lm(clineff ~ surreff, data = trialdat_fphr, weights = w)
r2label_fphr <- round(summary(trialmod_fphr)$r.squared, 2)

# Leave-one-out analysis
loo_preds_fphr <- rep(NA, 8)
for (i in 1:8) {

  trialmod_fphr <- lm(clineff ~ surreff,
                       data = trialdat_fphr[-i, ], weights = w)
  loo_preds_fphr[i] <- predict(trialmod_fphr,
                                newdata = trialdat_fphr[i, ])

}





# Make combined plot
# Make full data for figure
trialdat_rev$type <- "REVEAL 2.0"
trialdat_revlite$type <- "REVEAL Lite 2"
trialdat_comp$type <- "COMPERA"
trialdat_comp2$type <- "COMPERA 2.0"
trialdat_fphr$type <- "FPHR"
trialdat_full <- rbind(trialdat_rev, trialdat_revlite, trialdat_comp,
                       trialdat_comp2, trialdat_fphr)
trialdat_full$type <- factor(trialdat_full$type,
                             levels = c("REVEAL 2.0", "REVEAL Lite 2",
                                        "COMPERA", "COMPERA 2.0", "FPHR"))

r2labels <-
  data.frame(r2label = c(r2label_rev, r2label_revlite, r2label_comp,
                         r2label_comp2, r2label_fphr),
             type = factor(c("REVEAL 2.0", "REVEAL Lite 2", "COMPERA",
                      "COMPERA 2.0", "FPHR"),
                      levels = c("REVEAL 2.0", "REVEAL Lite 2", "COMPERA",
                                 "COMPERA 2.0", "FPHR")))

# Make figure, export as 7 x 5.25 pdf
ggplot(trialdat_full, aes(x = surreff, y = clineff)) +
  geom_point(aes(size = w), shape = 21) +
  facet_wrap(~type, nrow = 2) +
  xlab("Treatment effect on low risk status, log(OR)") +
  ylab("Treatment effect on clinical worsening, log(HR)") +
  ylim(-1, 0.1) +
  geom_smooth(method = 'lm', mapping = aes(weight = w), se = FALSE,
              color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  #annotate("text", x = 0.6, y = -0.75, label = r2label)
  geom_text(data = r2labels, size = 3.25,
            aes(x = 0.5, y = -0.875,
                label = paste("R^2", " == ", r2label, sep = "")),
            parse = TRUE) +
  theme_bw() +
  guides(size = "none")

# Make leave-one-out table
loo_tab <-
  data.frame("Observed trial log(HR)" = round(exp(trialdat_rev$clineff), 2),
             "REVEAL 2.0" = round(exp(loo_preds_rev), 2),
             "REVEAL Lite 2" = round(exp(loo_preds_revlite), 2),
             "COMPERA" = round(exp(loo_preds_comp), 2),
             "COMPERA 2.0" = round(exp(loo_preds_comp2), 2),
             "FPHR (non-inv)" = round(c(exp(loo_preds_fphr[1:2]), NA, NA, NA,
                                        exp(loo_preds_fphr[3:5])), 2))

# Output table
ft <- flextable(loo_tab) %>%
  theme_zebra()
#ft <- width(ft, j = 1, width = 2.5)

# Export table to Word
read_docx() %>%
  body_add_par("Meta-analysis leave-one-out table") %>%
  body_add_flextable(value = ft) %>%
  print(target = "meta_analysis_table.docx")








##########################################################################
# Sensitivity analysis using accelerated failure-time models

###########################################################################
# REVEAL 2.0 Meta-analyses
# Trial-region specific effects and associations
coef_rev_med <- list()
coef_rev_out <- list()
for (i in 1:8) {
  
  med_mod <- glm(reveal2_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- survreg(Surv(cw_day_full, cw_bin) ~ trt, dist = "weibull",
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  coef_rev_med[[i]] <- summary(med_mod)$coefficients
  coef_rev_out[[i]] <- summary(out_mod)$table
  
}

trialdat_rev <-
  data.frame(surreff = c(coef_rev_med[[1]][2, 1], coef_rev_med[[2]][2, 1],
                         coef_rev_med[[3]][2, 1], coef_rev_med[[4]][2, 1],
                         coef_rev_med[[5]][2, 1], coef_rev_med[[6]][2, 1],
                         coef_rev_med[[7]][2, 1], coef_rev_med[[8]][2, 1]),
             clineff = c(coef_rev_out[[1]][2, 1], coef_rev_out[[2]][2, 1],
                         coef_rev_out[[3]][2, 1], coef_rev_out[[4]][2, 1],
                         coef_rev_out[[5]][2, 1], coef_rev_out[[6]][2, 1],
                         coef_rev_out[[7]][2, 1], coef_rev_out[[8]][2, 1]),
             w = c(1/coef_rev_out[[1]][2, 2]^2, 1/coef_rev_out[[2]][2, 2]^2,
                   1/coef_rev_out[[3]][2, 2]^2, 1/coef_rev_out[[4]][2, 2]^2,
                   1/coef_rev_out[[5]][2, 2]^2, 1/coef_rev_out[[6]][2, 2]^2,
                   1/coef_rev_out[[7]][2, 2]^2, 1/coef_rev_out[[8]][2, 2]^2),
             n = c(table(dat$studyregion)))

trialmod_rev <- lm(clineff ~ surreff, data = trialdat_rev, weights = w)
r2label_rev <- round(summary(trialmod_rev)$r.squared, 2)

# Leave-one-out analysis
loo_preds_rev <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_rev <- lm(clineff ~ surreff, data = trialdat_rev[-i, ], weights = w)
  loo_preds_rev[i] <- predict(trialmod_rev, newdata = trialdat_rev[i, ])
  
}


###########################################################################
# REVEAL Lite Meta-analyses
# Trial-region specific effects and associations
coef_revlite_med <- list()
coef_revlite_out <- list()
for (i in 1:8) {
  
  med_mod <- glm(reveal_lite_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- survreg(Surv(cw_day_full, cw_bin) ~ trt, dist = "weibull",
                  data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_revlite_med[[i]] <- summary(med_mod)$coefficients
  coef_revlite_out[[i]] <- summary(out_mod)$table
  
}

trialdat_revlite <-
  data.frame(surreff = c(coef_revlite_med[[1]][2, 1],
                         coef_revlite_med[[2]][2, 1],
                         coef_revlite_med[[3]][2, 1],
                         coef_revlite_med[[4]][2, 1],
                         coef_revlite_med[[5]][2, 1],
                         coef_revlite_med[[6]][2, 1],
                         coef_revlite_med[[7]][2, 1],
                         coef_revlite_med[[8]][2, 1]),
             clineff = c(coef_revlite_out[[1]][2, 1],
                         coef_revlite_out[[2]][2, 1],
                         coef_revlite_out[[3]][2, 1],
                         coef_revlite_out[[4]][2, 1],
                         coef_revlite_out[[5]][2, 1],
                         coef_revlite_out[[6]][2, 1],
                         coef_revlite_out[[7]][2, 1],
                         coef_revlite_out[[8]][2, 1]),
             w = c(1/coef_revlite_out[[1]][2, 2]^2,
                   1/coef_revlite_out[[2]][2, 2]^2,
                   1/coef_revlite_out[[3]][2, 2]^2,
                   1/coef_revlite_out[[4]][2, 2]^2,
                   1/coef_revlite_out[[5]][2, 2]^2,
                   1/coef_revlite_out[[6]][2, 2]^2,
                   1/coef_revlite_out[[7]][2, 2]^2,
                   1/coef_revlite_out[[8]][2, 2]^2),
             n = c(table(dat$studyregion)))

trialmod_revlite <- lm(clineff ~ surreff, data = trialdat_revlite, weights = w)
r2label_revlite <- round(summary(trialmod_revlite)$r.squared, 2)

# Leave-one-out analysis
loo_preds_revlite <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_revlite <- lm(clineff ~ surreff,
                         data = trialdat_revlite[-i, ], weights = w)
  loo_preds_revlite[i] <- predict(trialmod_revlite,
                                  newdata = trialdat_revlite[i, ])
  
}


###########################################################################
# COMPERA Meta-analyses
# Trial-region specific effects and associations
coef_comp_med <- list()
coef_comp_out <- list()
for (i in 1:8) {
  
  med_mod <- glm(compera_full_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- survreg(Surv(cw_day_full, cw_bin) ~ trt, dist = "weibull",
                  data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_comp_med[[i]] <- summary(med_mod)$coefficients
  coef_comp_out[[i]] <- summary(out_mod)$table
  
}

trialdat_comp <-
  data.frame(surreff = c(coef_comp_med[[1]][2, 1], coef_comp_med[[2]][2, 1],
                         coef_comp_med[[3]][2, 1], coef_comp_med[[4]][2, 1],
                         coef_comp_med[[5]][2, 1], coef_comp_med[[6]][2, 1],
                         coef_comp_med[[7]][2, 1], coef_comp_med[[8]][2, 1]),
             clineff = c(coef_comp_out[[1]][2, 1], coef_comp_out[[2]][2, 1],
                         coef_comp_out[[3]][2, 1], coef_comp_out[[4]][2, 1],
                         coef_comp_out[[5]][2, 1], coef_comp_out[[6]][2, 1],
                         coef_comp_out[[7]][2, 1], coef_comp_out[[8]][2, 1]),
             w = c(1/coef_comp_out[[1]][2, 2]^2, 1/coef_comp_out[[2]][2, 2]^2,
                   1/coef_comp_out[[3]][2, 2]^2, 1/coef_comp_out[[4]][2, 2]^2,
                   1/coef_comp_out[[5]][2, 2]^2, 1/coef_comp_out[[6]][2, 2]^2,
                   1/coef_comp_out[[7]][2, 2]^2, 1/coef_comp_out[[8]][2, 2]^2),
             n = c(table(dat$studyregion)))

trialmod_comp <- lm(clineff ~ surreff, data = trialdat_comp, weights = w)
r2label_comp <- round(summary(trialmod_comp)$r.squared, 2)

# Leave-one-out analysis
loo_preds_comp <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_comp <- lm(clineff ~ surreff,
                      data = trialdat_comp[-i, ], weights = w)
  loo_preds_comp[i] <- predict(trialmod_comp,
                               newdata = trialdat_comp[i, ])
  
}


###########################################################################
# COMPERA 2.0 Meta-analyses
# Trial-region specific effects and associations
coef_comp2_med <- list()
coef_comp2_out <- list()
for (i in 1:8) {
  
  med_mod <- glm(compera_2_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- survreg(Surv(cw_day_full, cw_bin) ~ trt, dist = "weibull",
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_comp2_med[[i]] <- summary(med_mod)$coefficients
  coef_comp2_out[[i]] <- summary(out_mod)$table
  
}

trialdat_comp2 <-
  data.frame(surreff = c(coef_comp2_med[[1]][2, 1], coef_comp2_med[[2]][2, 1],
                         coef_comp2_med[[3]][2, 1], coef_comp2_med[[4]][2, 1],
                         coef_comp2_med[[5]][2, 1], coef_comp2_med[[6]][2, 1],
                         coef_comp2_med[[7]][2, 1], coef_comp2_med[[8]][2, 1]),
             clineff = c(coef_comp2_out[[1]][2, 1], coef_comp2_out[[2]][2, 1],
                         coef_comp2_out[[3]][2, 1], coef_comp2_out[[4]][2, 1],
                         coef_comp2_out[[5]][2, 1], coef_comp2_out[[6]][2, 1],
                         coef_comp2_out[[7]][2, 1], coef_comp2_out[[8]][2, 1]),
             w = c(1/coef_comp2_out[[1]][2, 2]^2,
                   1/coef_comp2_out[[2]][2, 2]^2,
                   1/coef_comp2_out[[3]][2, 2]^2,
                   1/coef_comp2_out[[4]][2, 2]^2,
                   1/coef_comp2_out[[5]][2, 2]^2,
                   1/coef_comp2_out[[6]][2, 2]^2,
                   1/coef_comp2_out[[7]][2, 2]^2,
                   1/coef_comp2_out[[8]][2, 2]^2),
             n = c(table(dat$studyregion)))

trialmod_comp2 <- lm(clineff ~ surreff, data = trialdat_comp2, weights = w)
r2label_comp2 <- round(summary(trialmod_comp2)$r.squared, 2)

# Leave-one-out analysis
loo_preds_comp2 <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_comp2 <- lm(clineff ~ surreff,
                       data = trialdat_comp2[-i, ], weights = w)
  loo_preds_comp2[i] <- predict(trialmod_comp2,
                                newdata = trialdat_comp2[i, ])
  
}


###########################################################################
# FPHR Meta-analyses
# Trial-region specific effects and associations
coef_fphr_med <- list()
coef_fphr_out <- list()
for (i in c(1, 2, 6, 7, 8)) {
  
  med_mod <- glm(fphr_noninv_low ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ],
                 family = "binomial")
  out_mod <- survreg(Surv(cw_day_full, cw_bin) ~ trt, dist = "weibull",
                  data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_fphr_med[[i]] <- summary(med_mod)$coefficients
  coef_fphr_out[[i]] <- summary(out_mod)$table
  
}

trialdat_fphr <-
  data.frame(surreff = c(coef_fphr_med[[1]][2, 1], coef_fphr_med[[2]][2, 1],
                         coef_fphr_med[[6]][2, 1], coef_fphr_med[[7]][2, 1],
                         coef_fphr_med[[8]][2, 1]),
             clineff = c(coef_fphr_out[[1]][2, 1], coef_fphr_out[[2]][2, 1],
                         coef_fphr_out[[6]][2, 1], coef_fphr_out[[7]][2, 1],
                         coef_fphr_out[[8]][2, 1]),
             w = c(1/coef_fphr_out[[1]][2, 2]^2,
                   1/coef_fphr_out[[2]][2, 2]^2,
                   1/coef_fphr_out[[6]][2, 2]^2,
                   1/coef_fphr_out[[7]][2, 2]^2,
                   1/coef_fphr_out[[8]][2, 2]^2),
             n = c(table(dat$studyregion)[c(1, 2, 6, 7, 8)]))

trialmod_fphr <- lm(clineff ~ surreff, data = trialdat_fphr, weights = w)
r2label_fphr <- round(summary(trialmod_fphr)$r.squared, 2)

# Leave-one-out analysis
loo_preds_fphr <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_fphr <- lm(clineff ~ surreff,
                      data = trialdat_fphr[-i, ], weights = w)
  loo_preds_fphr[i] <- predict(trialmod_fphr,
                               newdata = trialdat_fphr[i, ])
  
}





# Make combined plot
# Make full data for figure
trialdat_rev$type <- "REVEAL 2.0"
trialdat_revlite$type <- "REVEAL Lite 2"
trialdat_comp$type <- "COMPERA"
trialdat_comp2$type <- "COMPERA 2.0"
trialdat_fphr$type <- "FPHR"
trialdat_full <- rbind(trialdat_rev, trialdat_revlite, trialdat_comp,
                       trialdat_comp2, trialdat_fphr)
trialdat_full$type <- factor(trialdat_full$type,
                             levels = c("REVEAL 2.0", "REVEAL Lite 2",
                                        "COMPERA", "COMPERA 2.0", "FPHR"))

r2labels <-
  data.frame(r2label = c(r2label_rev, r2label_revlite, r2label_comp,
                         r2label_comp2, r2label_fphr),
             type = factor(c("REVEAL 2.0", "REVEAL Lite 2", "COMPERA",
                             "COMPERA 2.0", "FPHR"),
                           levels = c("REVEAL 2.0", "REVEAL Lite 2",
                                      "COMPERA", "COMPERA 2.0", "FPHR")))

# Make figure
ggplot(trialdat_full, aes(x = surreff, y = clineff)) +
  geom_point(aes(size = w), shape = 21) +
  facet_wrap(~type, nrow = 2) +
  xlab("Treatment effect on low risk status, log(OR)") +
  ylab("Treatment effect on clinical worsening") +
  ylim(-0.1, 0.9) +
  geom_smooth(method = 'lm', mapping = aes(weight = w), se = FALSE,
              color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  #annotate("text", x = 0.6, y = -0.75, label = r2label)
  geom_text(data = r2labels, size = 3.25,
            aes(x = 0.5, y = 0.85,
                label = paste("R^2", " == ", r2label, sep = "")),
            parse = TRUE) +
  theme_bw() +
  guides(size = "none")


# Make leave-one-out table
loo_tab <-
  data.frame("Observed trial log(HR)" = round(trialdat_rev$clineff, 2),
             "REVEAL 2.0" = round(loo_preds_rev, 2),
             "REVEAL Lite 2" = round(loo_preds_revlite, 2),
             "COMPERA" = round(loo_preds_comp, 2),
             "COMPERA 2.0" = round(loo_preds_comp2, 2),
             "FPHR (non-inv)" = round(c(loo_preds_fphr[1:2], NA, NA, NA,
                                        loo_preds_fphr[3:5]), 2))

# Output table
ft2 <- flextable(loo_tab) %>%
  theme_zebra()
#ft2 <- width(ft, j = 1, width = 2.5)

# Export table to Word
read_docx() %>%
  body_add_par("Meta-analysis leave-one-out table") %>%
  body_add_flextable(value = ft2) %>%
  print(target = "meta_analysis_table2.docx")






############################################################################
# Sensitivity analysis: Use ordinal risk score values

###########################################################################
# REVEAL 2.0 Meta-analyses
# Trial-region specific effects and associations
coef_rev_med <- list()
coef_rev_out <- list()
for (i in 1:8) {
  
  med_mod <- glm(reveal2_wk16 ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_rev_med[[i]] <- summary(med_mod)$coefficients
  coef_rev_out[[i]] <- summary(out_mod)$coefficients
  
}

trialdat_rev <-
  data.frame(surreff = c(coef_rev_med[[1]][2, 1], coef_rev_med[[2]][2, 1],
                         coef_rev_med[[3]][2, 1], coef_rev_med[[4]][2, 1],
                         coef_rev_med[[5]][2, 1], coef_rev_med[[6]][2, 1],
                         coef_rev_med[[7]][2, 1], coef_rev_med[[8]][2, 1]),
             clineff = c(coef_rev_out[[1]][1, 1], coef_rev_out[[2]][1, 1],
                         coef_rev_out[[3]][1, 1], coef_rev_out[[4]][1, 1],
                         coef_rev_out[[5]][1, 1], coef_rev_out[[6]][1, 1],
                         coef_rev_out[[7]][1, 1], coef_rev_out[[8]][1, 1]),
             w = c(1/coef_rev_out[[1]][1, 3]^2, 1/coef_rev_out[[2]][1, 3]^2,
                   1/coef_rev_out[[3]][1, 3]^2, 1/coef_rev_out[[4]][1, 3]^2,
                   1/coef_rev_out[[5]][1, 3]^2, 1/coef_rev_out[[6]][1, 3]^2,
                   1/coef_rev_out[[7]][1, 3]^2, 1/coef_rev_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_rev <- lm(clineff ~ surreff, data = trialdat_rev, weights = w)
r2label_rev <- round(summary(trialmod_rev)$r.squared, 2)

# Leave-one-out analysis
loo_preds_rev <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_rev <- lm(clineff ~ surreff, data = trialdat_rev[-i, ], weights = w)
  loo_preds_rev[i] <- predict(trialmod_rev, newdata = trialdat_rev[i, ])
  
}


###########################################################################
# REVEAL Lite Meta-analyses
# Trial-region specific effects and associations
coef_revlite_med <- list()
coef_revlite_out <- list()
for (i in 1:8) {
  
  med_mod <- glm(reveal_lite_wk16 ~ trt,
                 data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_revlite_med[[i]] <- summary(med_mod)$coefficients
  coef_revlite_out[[i]] <- summary(out_mod)$coefficients
  
}

trialdat_revlite <-
  data.frame(surreff = c(coef_revlite_med[[1]][2, 1],
                         coef_revlite_med[[2]][2, 1],
                         coef_revlite_med[[3]][2, 1],
                         coef_revlite_med[[4]][2, 1],
                         coef_revlite_med[[5]][2, 1],
                         coef_revlite_med[[6]][2, 1],
                         coef_revlite_med[[7]][2, 1],
                         coef_revlite_med[[8]][2, 1]),
             clineff = c(coef_revlite_out[[1]][1, 1],
                         coef_revlite_out[[2]][1, 1],
                         coef_revlite_out[[3]][1, 1],
                         coef_revlite_out[[4]][1, 1],
                         coef_revlite_out[[5]][1, 1],
                         coef_revlite_out[[6]][1, 1],
                         coef_revlite_out[[7]][1, 1],
                         coef_revlite_out[[8]][1, 1]),
             w = c(1/coef_revlite_out[[1]][1, 3]^2,
                   1/coef_revlite_out[[2]][1, 3]^2,
                   1/coef_revlite_out[[3]][1, 3]^2,
                   1/coef_revlite_out[[4]][1, 3]^2,
                   1/coef_revlite_out[[5]][1, 3]^2,
                   1/coef_revlite_out[[6]][1, 3]^2,
                   1/coef_revlite_out[[7]][1, 3]^2,
                   1/coef_revlite_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_revlite <- lm(clineff ~ surreff, data = trialdat_revlite, weights = w)
r2label_revlite <- round(summary(trialmod_revlite)$r.squared, 2)

# Leave-one-out analysis
loo_preds_revlite <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_revlite <- lm(clineff ~ surreff,
                         data = trialdat_revlite[-i, ], weights = w)
  loo_preds_revlite[i] <- predict(trialmod_revlite,
                                  newdata = trialdat_revlite[i, ])
  
}


###########################################################################
# COMPERA Meta-analyses
# Trial-region specific effects and associations
dat$compera_cat0 <- round(dat$compera_full_wk0 + 0.01)
dat$compera_cat16 <- round(dat$compera_full_wk16 + 0.01)
coef_comp_med <- list()
coef_comp_out <- list()
for (i in 1:8) {
  
  med_mod <- polr(factor(compera_cat16) ~ trt,
                  data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_comp_med[[i]] <- summary(med_mod)$coefficients
  coef_comp_out[[i]] <- summary(out_mod)$coefficients
  
}

trialdat_comp <-
  data.frame(surreff = c(coef_comp_med[[1]][1, 1], coef_comp_med[[2]][1, 1],
                         coef_comp_med[[3]][1, 1], coef_comp_med[[4]][1, 1],
                         coef_comp_med[[5]][1, 1], coef_comp_med[[6]][1, 1],
                         coef_comp_med[[7]][1, 1], coef_comp_med[[8]][1, 1]),
             clineff = c(coef_comp_out[[1]][1, 1], coef_comp_out[[2]][1, 1],
                         coef_comp_out[[3]][1, 1], coef_comp_out[[4]][1, 1],
                         coef_comp_out[[5]][1, 1], coef_comp_out[[6]][1, 1],
                         coef_comp_out[[7]][1, 1], coef_comp_out[[8]][1, 1]),
             w = c(1/coef_comp_out[[1]][1, 3]^2, 1/coef_comp_out[[2]][1, 3]^2,
                   1/coef_comp_out[[3]][1, 3]^2, 1/coef_comp_out[[4]][1, 3]^2,
                   1/coef_comp_out[[5]][1, 3]^2, 1/coef_comp_out[[6]][1, 3]^2,
                   1/coef_comp_out[[7]][1, 3]^2, 1/coef_comp_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_comp <- lm(clineff ~ surreff, data = trialdat_comp, weights = w)
r2label_comp <- round(summary(trialmod_comp)$r.squared, 2)

# Leave-one-out analysis
loo_preds_comp <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_comp <- lm(clineff ~ surreff,
                      data = trialdat_comp[-i, ], weights = w)
  loo_preds_comp[i] <- predict(trialmod_comp,
                               newdata = trialdat_comp[i, ])
  
}


###########################################################################
# COMPERA 2.0 Meta-analyses
# Trial-region specific effects and associations
dat$compera2_cat0 <- round(dat$compera_2_wk0 + 0.01)
dat$compera2_cat16 <- round(dat$compera_2_wk16 + 0.01)
coef_comp2_med <- list()
coef_comp2_out <- list()
for (i in 1:8) {
  
  med_mod <- polr(factor(compera2_cat16) ~ trt,
                  data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_comp2_med[[i]] <- summary(med_mod)$coefficients
  coef_comp2_out[[i]] <- summary(out_mod)$coefficients
  
}

trialdat_comp2 <-
  data.frame(surreff = c(coef_comp2_med[[1]][1, 1], coef_comp2_med[[2]][1, 1],
                         coef_comp2_med[[3]][1, 1], coef_comp2_med[[4]][1, 1],
                         coef_comp2_med[[5]][1, 1], coef_comp2_med[[6]][1, 1],
                         coef_comp2_med[[7]][1, 1], coef_comp2_med[[8]][1, 1]),
             clineff = c(coef_comp2_out[[1]][1, 1], coef_comp2_out[[2]][1, 1],
                         coef_comp2_out[[3]][1, 1], coef_comp2_out[[4]][1, 1],
                         coef_comp2_out[[5]][1, 1], coef_comp2_out[[6]][1, 1],
                         coef_comp2_out[[7]][1, 1], coef_comp2_out[[8]][1, 1]),
             w = c(1/coef_comp2_out[[1]][1, 3]^2,
                   1/coef_comp2_out[[2]][1, 3]^2,
                   1/coef_comp2_out[[3]][1, 3]^2,
                   1/coef_comp2_out[[4]][1, 3]^2,
                   1/coef_comp2_out[[5]][1, 3]^2,
                   1/coef_comp2_out[[6]][1, 3]^2,
                   1/coef_comp2_out[[7]][1, 3]^2,
                   1/coef_comp2_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)))

trialmod_comp2 <- lm(clineff ~ surreff, data = trialdat_comp2, weights = w)
r2label_comp2 <- round(summary(trialmod_comp2)$r.squared, 2)

# Leave-one-out analysis
loo_preds_comp2 <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_comp2 <- lm(clineff ~ surreff,
                       data = trialdat_comp2[-i, ], weights = w)
  loo_preds_comp2[i] <- predict(trialmod_comp2,
                                newdata = trialdat_comp2[i, ])
  
}


###########################################################################
# FPHR Meta-analyses
# Trial-region specific effects and associations
dat$fphr_cat16 <- 3 - dat$fphr_noninv_wk16
dat$fphr_cat16[dat$fphr_cat16 == 0] <- 1
dat$fphr_cat16 <- factor(dat$fphr_cat16)
levels(dat$fphr_cat16) <- c("Low", "Intermediate", "High")
coef_fphr_med <- list()
coef_fphr_out <- list()
for (i in c(1, 2, 6, 7, 8)) {
  
  med_mod <- polr(factor(fphr_cat16) ~ trt,
                  data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  out_mod <- coxph(Surv(cw_day_full, cw_bin) ~ trt,
                   data = dat[dat$studyregion == unique(dat$studyregion)[i], ])
  
  coef_fphr_med[[i]] <- summary(med_mod)$coefficients
  coef_fphr_out[[i]] <- summary(out_mod)$coefficients
  
}

trialdat_fphr <-
  data.frame(surreff = c(coef_fphr_med[[1]][1, 1], coef_fphr_med[[2]][1, 1],
                         coef_fphr_med[[6]][1, 1], coef_fphr_med[[7]][1, 1],
                         coef_fphr_med[[8]][1, 1]),
             clineff = c(coef_fphr_out[[1]][1, 1], coef_fphr_out[[2]][1, 1],
                         coef_fphr_out[[6]][1, 1], coef_fphr_out[[7]][1, 1],
                         coef_fphr_out[[8]][1, 1]),
             w = c(1/coef_fphr_out[[1]][1, 3]^2,
                   1/coef_fphr_out[[2]][1, 3]^2,
                   1/coef_fphr_out[[6]][1, 3]^2,
                   1/coef_fphr_out[[7]][1, 3]^2,
                   1/coef_fphr_out[[8]][1, 3]^2),
             n = c(table(dat$studyregion)[c(1, 2, 6, 7, 8)]))

trialmod_fphr <- lm(clineff ~ surreff, data = trialdat_fphr, weights = w)
r2label_fphr <- round(summary(trialmod_fphr)$r.squared, 2)

# Leave-one-out analysis
loo_preds_fphr <- rep(NA, 8)
for (i in 1:8) {
  
  trialmod_fphr <- lm(clineff ~ surreff,
                      data = trialdat_fphr[-i, ], weights = w)
  loo_preds_fphr[i] <- predict(trialmod_fphr,
                               newdata = trialdat_fphr[i, ])
  
}





# Make combined plot
# Make full data for figure
trialdat_rev$type <- "REVEAL 2.0"
trialdat_revlite$type <- "REVEAL Lite 2"
trialdat_comp$type <- "COMPERA"
trialdat_comp2$type <- "COMPERA 2.0"
trialdat_fphr$type <- "FPHR"
trialdat_full <- rbind(trialdat_rev, trialdat_revlite, trialdat_comp,
                       trialdat_comp2, trialdat_fphr)
trialdat_full$type <- factor(trialdat_full$type,
                             levels = c("REVEAL 2.0", "REVEAL Lite 2",
                                        "COMPERA", "COMPERA 2.0", "FPHR"))

r2labels <-
  data.frame(r2label = c(r2label_rev, r2label_revlite, r2label_comp,
                         r2label_comp2, r2label_fphr),
             type = factor(c("REVEAL 2.0", "REVEAL Lite 2", "COMPERA",
                             "COMPERA 2.0", "FPHR"),
                           levels = c("REVEAL 2.0", "REVEAL Lite 2",
                                      "COMPERA", "COMPERA 2.0", "FPHR")))

# Make figure
ggplot(trialdat_full, aes(x = surreff, y = clineff)) +
  geom_point(aes(size = w), shape = 21) +
  facet_wrap(~type, nrow = 2) +
  xlab("Treatment effect on ordinal risk score, log(OR)") +
  ylab("Treatment effect on clinical worsening, log(HR)") +
  ylim(-1, 0.1) +
  geom_smooth(method = 'lm', mapping = aes(weight = w), se = FALSE,
              color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  #annotate("text", x = 0.6, y = -0.75, label = r2label)
  geom_text(data = r2labels, size = 3.25,
            aes(x = 0, y = -0.875,
                label = paste("R^2", " == ", r2label, sep = "")),
            parse = TRUE) +
  theme_bw() +
  guides(size = "none")

# Make leave-one-out table
loo_tab <-
  data.frame("Observed trial log(HR)" = round(exp(trialdat_rev$clineff), 2),
             "REVEAL 2.0" = round(exp(loo_preds_rev), 2),
             "REVEAL Lite 2" = round(exp(loo_preds_revlite), 2),
             "COMPERA" = round(exp(loo_preds_comp), 2),
             "COMPERA 2.0" = round(exp(loo_preds_comp2), 2),
             "FPHR (non-inv)" = round(c(exp(loo_preds_fphr[1:2]), NA, NA, NA,
                                        exp(loo_preds_fphr[3:5])), 2))

# Output table
ft <- flextable(loo_tab) %>%
  theme_zebra()
#ft <- width(ft, j = 1, width = 2.5)

# Export table to Word
read_docx() %>%
  body_add_par("Meta-analysis leave-one-out table") %>%
  body_add_flextable(value = ft) %>%
  print(target = "meta_analysis_table_old.docx")