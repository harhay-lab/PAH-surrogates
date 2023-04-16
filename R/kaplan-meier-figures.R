# Clear data and load libraries
rm(list = ls())
library(ggsurvfit)
library(glue)
library(gridExtra)
library(cowplot)

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
dat$trt[dat$trt == "Active"] <- "Experimental"

# Make low risk variables
# Low REVEAL 2.0
dat$reveal_cat16 <- NA
dat$reveal_cat16[round(dat$reveal2_wk16) <= 6] <- 1
dat$reveal_cat16[round(dat$reveal2_wk16) >= 7] <- 0

# Categorical REVEAL Lite
dat$reveal_lite_cat16 <- NA
dat$reveal_lite_cat16[round(dat$reveal_lite_wk16) <= 5] <- 1
dat$reveal_lite_cat16[round(dat$reveal_lite_wk16) >= 6] <- 0

# Categorical COMPERA (round to avoid R issues)
dat$compera_cat16 <- NA
dat$compera_cat16[round(dat$compera_full_wk16 + 0.01) == 1] <- 1
dat$compera_cat16[round(dat$compera_full_wk16 + 0.01) != 1] <- 0

# Categorical COMPERA 2.0
dat$compera2_cat16 <- NA
dat$compera2_cat16[round(dat$compera_2_wk16 + 0.01) == 1] <- 1
dat$compera2_cat16[round(dat$compera_2_wk16 + 0.01) != 1] <- 0

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
p1 <- survfit2(Surv(cw_day_full2, cw_bin) ~ trt,
               data = dat[dat$cw_day_full > 16*7, ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "A", x = "Weeks",
       y = "Cumulative survival free from clinical worsening") +
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250, y = 0.67, size = 3,
    label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full2,
                       cw_bin) ~ trt,
                       data = dat[dat$cw_day_full > 16*7, ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(breaks = c("Control", "Experimental"),
                      values = c("grey", "black"), name = "Treatment arm") +
  scale_fill_manual(breaks = c("Control", "Experimental"),
                    values = c("grey", "black"), name = "Treatment arm") +
  scale_linetype_manual(values = 1:2,
                        labels = c("Control", "Experimental"),
                        name = "Treatment arm") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.25, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 9))

p2 <- survfit2(Surv(death_day_full2, death_bin) ~ trt,
               data = dat[dat$death_day_full > 16*7, ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "B", x = "Weeks",
       y = "Cumulative survival") + 
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250, y = 0.75, size = 3,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full2,
                              death_bin) ~ trt,
                              data = dat[dat$death_day_full > 16*7, ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(breaks = c("Control", "Experimental"),
                      values = c("grey", "black"), name = "Treatment arm") +
  scale_fill_manual(breaks = c("Control", "Experimental"),
                    values = c("grey", "black"), name = "Treatment arm") +
  scale_linetype_manual(values = 1:2,
                        labels = c("Control", "Experimental"),
                        name = "Treatment arm") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.25, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 9))

# COMPERA panels
p3 <- survfit2(Surv(cw_day_full2, cw_bin) ~ compera_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$compera_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "C", x = "Weeks",
       y = "Cumulative survival free from clinical worsening") + 
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250-16*7, y = 0.85, size = 3,
           label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full2,
                              cw_bin) ~ compera_cat16,
                              data = dat[dat$cw_day_full > 16*7 &
                                         !is.na(dat$compera_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "COMPERA score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "COMPERA score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "COMPERA score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.25, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 9))

p4 <- survfit2(Surv(death_day_full2, death_bin) ~ compera_cat16,
               data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$compera_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "D", x = "Weeks",
       y = "Cumulative survival") + 
  ylim(c(0.25, 1)) +
  annotate("text", x = 1250-16*7, y = 0.7, size = 3,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full2,
                              death_bin) ~ compera_cat16,
                              data = dat[dat$death_day_full > 16*7 &
                                         !is.na(dat$compera_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "COMPERA score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "COMPERA score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "COMPERA score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.25, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 9))

# Make figure for paper (edit dimensions as needed)
pdf("Output/kaplan-meier-figure.pdf", height = 7.5, width = 6)
#grid.arrange(p1, p2, p3, p4, nrow = 2) # doesn't add risk tables
plot_grid(ggsurvfit_build(p1), ggsurvfit_build(p2),
          ggsurvfit_build(p3), ggsurvfit_build(p4), ncol = 2)
dev.off()






###########################################################################
# Figure for supplement (time to CW/survival after other risk scores)

# REVEAL 2.0 panels
p5 <- survfit2(Surv(cw_day_full2, cw_bin) ~ reveal_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$reveal_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "A", x = "Weeks",
       y = "Cumulative survival free from CW") + 
  ylim(c(0.2, 1)) +
  annotate("text", x = 1250-16*7, y = 0.85, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full2,
                              cw_bin) ~ reveal_cat16,
                              data = dat[dat$cw_day_full > 16*7 &
                                         !is.na(dat$reveal_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "REVEAL 2.0 score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "REVEAL 2.0 score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "REVEAL 2.0 score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

p6 <- survfit2(Surv(death_day_full2, death_bin) ~ reveal_cat16,
               data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$reveal_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "A", x = "Weeks",
       y = "Cumulative survival") + 
  ylim(c(0.2, 1)) +
  annotate("text", x = 1250-16*7, y = 0.7, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full2,
                              death_bin) ~ reveal_cat16,
                              data = dat[dat$death_day_full > 16*7 &
                                         !is.na(dat$reveal_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "REVEAL 2.0 score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "REVEAL 2.0 score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "REVEAL 2.0 score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

# REVEAL Lite panels
p7 <- survfit2(Surv(cw_day_full2, cw_bin) ~ reveal_lite_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$reveal_lite_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "B", x = "Weeks",
       y = "Cumulative survival free from CW") + 
  ylim(c(0.2, 1)) +
  annotate("text", x = 1250-16*7, y = 0.85, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full2,
                              cw_bin) ~ reveal_lite_cat16,
                              data = dat[dat$cw_day_full > 16*7 &
                                    !is.na(dat$reveal_lite_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "REVEAL Lite 2 score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "REVEAL Lite 2 score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "REVEAL Lite 2 score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

p8 <- survfit2(Surv(death_day_full2, death_bin) ~ reveal_lite_cat16,
               data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$reveal_lite_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "B", x = "Weeks",
       y = "Cumulative survival") + 
  ylim(c(0.2, 1)) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  annotate("text", x = 1250-16*7, y = 0.67, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full2,
                              death_bin) ~ reveal_lite_cat16,
                              data = dat[dat$death_day_full > 16*7 &
                                    !is.na(dat$reveal_lite_cat16), ]))}")) +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "REVEAL Lite 2 score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "REVEAL Lite 2 score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "REVEAL Lite 2 score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

# COMPERA 2.0 Panels
p9 <- survfit2(Surv(cw_day_full2, cw_bin) ~ compera2_cat16,
               data = dat[dat$cw_day_full > 16*7 &
                            !is.na(dat$compera2_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "C", x = "Weeks",
       y = "Cumulative survival free from CW") + 
  ylim(c(0.2, 1)) +
  annotate("text", x = 1250-16*7, y = 0.88, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full2,
                              cw_bin) ~ compera2_cat16,
                              data = dat[dat$cw_day_full > 16*7 &
                                         !is.na(dat$compera2_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "COMPERA 2.0 score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "COMPERA 2.0 score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "COMPERA 2.0 score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

p10 <- survfit2(Surv(death_day_full2, death_bin) ~ compera2_cat16,
                data = dat[dat$death_day_full > 16*7 &
                            !is.na(dat$compera2_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "C", x = "Weeks",
       y = "Cumulative survival") + 
  ylim(c(0.2, 1)) +
  annotate("text", x = 1250-16*7, y = 0.7, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full2,
                              death_bin) ~ compera2_cat16,
                              data = dat[dat$death_day_full > 16*7 &
                                         !is.na(dat$compera2_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "COMPERA 2.0 score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "COMPERA 2.0 score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "COMPERA 2.0 score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

# FPHR panels
p11 <- survfit2(Surv(cw_day_full2, cw_bin) ~ fphr_cat16,
                data = dat[dat$cw_day_full > 16*7 &
                             !is.na(dat$fphr_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "D", x = "Weeks",
       y = "Cumulative survival free from CW") +
  ylim(c(0.2, 1)) +
  annotate("text", x = 1250-16*7, y = 0.85, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(cw_day_full2,
                              cw_bin) ~ fphr_cat16,
                              data = dat[dat$cw_day_full > 16*7 &
                                         !is.na(dat$fphr_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "FPHR score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "FPHR score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "FPHR score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

p12 <- survfit2(Surv(death_day_full2, death_bin) ~ fphr_cat16,
                data = dat[dat$death_day_full > 16*7 &
                             !is.na(dat$fphr_cat16), ]) %>% 
  ggsurvfit(linetype_aes = TRUE) +
  labs(title = "D", x = "Weeks",
       y = "Cumulative survival") + 
  ylim(c(0.2, 1)) +
  annotate("text", x = 1250-16*7, y = 0.7, size = 9/.pt,
           label = glue::glue("{survfit2_p(survfit2(Surv(death_day_full2,
                              death_bin) ~ fphr_cat16,
                              data = dat[dat$death_day_full > 16*7 &
                                         !is.na(dat$fphr_cat16), ]))}")) +
  add_confidence_interval() +
  add_risktable(risktable_stats = c("n.risk")) +
  add_risktable_strata_symbol() +
  scale_colour_manual(values = c("grey", "black"),
                      labels = c("Not low risk", "Low risk"),
                      name = "FPHR score") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("Not low risk", "Low risk"),
                    name = "FPHR score") +
  scale_linetype_manual(values = 2:1,
                        labels = c("Not low risk", "Low risk"),
                        name = "FPHR score") +
  scale_x_continuous(breaks = c(0, (52-16)*7, (104-16)*7, (156-16)*7,
                                (208-16)*7),
                     labels = c(16, 52, 104, 156, 208)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_text(hjust = 0.5, size = 11),
        axis.title.y = element_text(hjust = 0.5, size = 11),
        legend.position = c(0.22, 0.27),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))


# Make supplement figures
pdf("Output/kaplan-meier-figure2.pdf", height = 7.5, width = 6)
#grid.arrange(p5, p7, p9, p11, nrow = 2)
plot_grid(ggsurvfit_build(p5), ggsurvfit_build(p7),
          ggsurvfit_build(p9), ggsurvfit_build(p11), ncol = 2)
dev.off()

pdf("Output/kaplan-meier-figure3.pdf", height = 7.5, width = 6)
#grid.arrange(p6, p8, p10, p12, nrow = 2)
plot_grid(ggsurvfit_build(p6), ggsurvfit_build(p8),
          ggsurvfit_build(p10), ggsurvfit_build(p12), ncol = 2)
dev.off()
