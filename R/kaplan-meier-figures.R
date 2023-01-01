# Clear data and load libraries
rm(list = ls())
library(ggsurvfit)
library(glue)
library(gridExtra)

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
dat$trt <- dat$study_arm_assignment

# Make low risk variables
# Low REVEAL 2.0
dat$reveal_cat16 <- NA
dat$reveal_cat16[round(dat$reveal2_wk16) <= 6] <- 1
dat$reveal_cat16[round(dat$reveal2_wk16) >= 7] <- 0

# Categorical REVEAL Lite
dat$reveal_lite_cat16 <- NA
dat$reveal_lite_cat16[round(dat$reveal_lite_wk16) <= 5] <- 1
dat$reveal_lite_cat16[round(dat$reveal_lite_wk16) >= 6] <- 0

# Categorical COMPERA
dat$compera_cat16 <- NA
dat$compera_cat16[round(dat$compera_full_wk16) == 1] <- 1
dat$compera_cat16[round(dat$compera_full_wk16) != 1] <- 0

# Categorical COMPERA 2.0
dat$compera2_cat16 <- NA
dat$compera2_cat16[round(dat$compera_2_wk16) == 1] <- 1
dat$compera2_cat16[round(dat$compera_2_wk16) != 1] <- 0

# Categorical FPHR (non-inv)
dat$fphr_cat16 <- NA
dat$fphr_cat16[dat$fphr_noninv_wk16 == 2 | dat$fphr_noninv_wk16 == 3] <- 1
dat$fphr_cat16[dat$fphr_noninv_wk16 == 0 | dat$fphr_noninv_wk16 == 1] <- 0

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

# Create variables for time since surrogate meaasurement
dat$cw_day_full2 <- dat$cw_day_full - 16*7
dat$death_day_full2 <- dat$death_day_full - 16*7

# Make figure panels
# Treatment panels
p1 <- survfit2(Surv(cw_day_full, cw_bin) ~ trt,
               data = dat) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "A", x = "Years",
       y = "Cumulative survival free from clinical worsening") +
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250, y = 0.67,
    label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full,
                       cw_bin) ~ trt, data = dat))}")) +
  add_confidence_interval() +
  scale_colour_manual(breaks = c("Active", "Control"),
                      values = c("grey", "black")) +
  scale_fill_manual(breaks = c("Active", "Control"),
                      values = c("grey", "black")) +
  scale_x_continuous(breaks = c(0, 365, 730, 1095, 1460), labels = 0:4) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(fill = "white", color = "black"))

p2 <- survfit2(Surv(death_day_full, death_bin) ~ trt,
               data = dat) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "B", x = "Years",
       y = "Cumulative survival") + 
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250, y = 0.75,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full,
                              death_bin) ~ trt, data = dat))}")) +
  add_confidence_interval() +
  scale_colour_manual(breaks = c("Active", "Control"),
                      values = c("grey", "black")) +
  scale_fill_manual(breaks = c("Active", "Control"),
                    values = c("grey", "black")) +
  scale_x_continuous(breaks = c(0, 365, 730, 1095, 1460), labels = 0:4) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(fill = "white", color = "black"))

# COMPERA panels
p3 <- survfit2(Surv(cw_day_full2, cw_bin) ~ compera_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$compera_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "C", x = "Years",
       y = "Cumulative survival free from clinical worsening") + 
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250-16*7, y = 0.85,
           label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full2,
                              cw_bin) ~ compera_cat16,
                              data = dat[dat$cw_day_full > 16*7 &
                                         !is.na(dat$compera_cat16), ]))}")) +
  add_confidence_interval() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk")) +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk")) +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk")) +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(fill = "white", color = "black"))

p4 <- survfit2(Surv(death_day_full2, death_bin) ~ compera_cat16,
               data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$compera_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "D", x = "Years",
       y = "Cumulative survival") + 
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250-16*7, y = 0.7,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full2,
                              death_bin) ~ compera_cat16,
                              data = dat[dat$death_day_full > 16*7 &
                                         !is.na(dat$compera_cat16), ]))}")) +
  add_confidence_interval() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk")) +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk")) +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk")) +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(fill = "white", color = "black"))

# Make figure for paper
pdf("Output/kaplan-meier-figure.pdf", height = 7.5, width = 6)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()






###########################################################################
# Figure for supplement (time to CW/survival after other risk scores)

# REVEAL 2.0 panels
p5 <- survfit2(Surv(cw_day_full2, cw_bin) ~ reveal_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$reveal_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "REVEAL 2.0", x = "Years",
       y = "Prop. survived w/o CW") + 
  add_confidence_interval() +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

p6 <- survfit2(Surv(death_day_full2, death_bin) ~ reveal_cat16,
               data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$reveal_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "REVEAL 2.0", x = "Years",
       y = "Prop. survived") + 
  add_confidence_interval() +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

# REVEAL Lite panels
p7 <- survfit2(Surv(cw_day_full2, cw_bin) ~ reveal_lite_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$reveal_lite_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "REVEAL Lite", x = "Years",
       y = "Prop. survived w/o CW") + 
  add_confidence_interval() +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

p8 <- survfit2(Surv(death_day_full2, death_bin) ~ reveal_lite_cat16,
               data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$reveal_lite_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "REVEAL Lite", x = "Years",
       y = "Prop. survived") + 
  add_confidence_interval() +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

# COMPERA 2.0 Panels
p9 <- survfit2(Surv(cw_day_full2, cw_bin) ~ compera2_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$compera2_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "COMPERA 2.0", x = "Years",
       y = "Prop. survived w/o CW") + 
  add_confidence_interval() +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

p10 <- survfit2(Surv(death_day_full2, death_bin) ~ compera2_cat16,
                data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$compera2_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "COMPERA 2.0", x = "Years",
       y = "Prop. survived") + 
  add_confidence_interval() +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

# FPHR panels
p11 <- survfit2(Surv(cw_day_full2, cw_bin) ~ fphr_cat16,
                data = dat[dat$cw_day_full > 16*7 &
                             !is.na(dat$fphr_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "FPHR", x = "Years",
       y = "Prop. survived w/o CW") + 
  add_confidence_interval() +
  scale_colour_manual(values = c("red", "blue"),
                      labels = c("Not low risk", "Low risk")) +
  scale_fill_manual(values = c("red", "blue"),
                      labels = c("Not low risk", "Low risk")) +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk")) +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  theme()

p12 <- survfit2(Surv(death_day_full2, death_bin) ~ fphr_cat16,
                data = dat[dat$death_day_full > 16*7 &
                             !is.na(dat$fphr_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "FPHR", x = "Years",
       y = "Prop. survived") + 
  add_confidence_interval() +
  scale_colour_manual(values = c("red", "blue"),
                      labels = c("Not low risk", "Low risk")) +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Not low risk", "Low risk")) +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk")) +
  scale_x_continuous(breaks = c(0, 365-16*7, 730-16*7, 1095-16*7, 1460-16*7),
                     labels = c("16 weeks", "1", "2", "3", "4")) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")


# Make supplement figure
pdf("kaplan-meier-figure2.pdf", height = 8.5, width = 6)
grid.arrange(p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)
dev.off()
