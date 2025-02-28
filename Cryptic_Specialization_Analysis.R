
# Cryptic specialization 
# Puerto Rico Transplants
# T. Lindsay PhD Chapter 2

# Set up  -----------------------------------------------------------------

# Environment 
set.seed(600)
setwd('~/Desktop/GITHUB/TL_Trans/')

# Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library("lubridate")
library(stringr)
library(hms)
library("survival")
library("survminer")
library(gridExtra)    # for grid.arrange
library(cowplot)      # arranging plots plot_grid()
library(grid)
library(ggfortify)    # pca plots 
library(flexsurv)     # for parametric survival curves 
library(ggpattern)    # for patterns
library(scales)       # scale_y_continuous 

# Data setup  -------------------------------------------------------------

# raw data files 
raw <- read.csv('TL_Trans_Results.csv')
PAM <- read.csv('PAM/TL_TRANS_PAM_Aug04_13.csv')
surv <- read.csv('Survival/TL_Trans_Survival_KP_Format.csv')
temp_DO <- read.csv('Meta_Environmental_Data/TLPR21_Temp_DO_Raw.csv') 
light_shallow <- read.csv('Meta_Environmental_Data/LIGHT_SHALLOW_ALL.csv') #%>% select(!c(SENSOR, old_date))
light_deep <- read.csv('Meta_Environmental_Data/LIGHT_DEEP_ALL.csv') #%>% select(!c(SENSOR, old_date))

# find mean of raw by genotype 
raw <- raw %>%
  group_by(colony_id, species, group, genotype, treatment, full_treatment) %>%
  summarise(A = mean(A),
            CA = mean(CA),
            di= mean(di),
            Cdi= mean(Cdi),
            D2 = mean(D2),
            Dsmall = mean(Dsmall),
            chla.ug.cm2= mean(chla.ug.cm2),
            chlc2.ug.cm2= mean(chlc2.ug.cm2),
            chla.sym= mean(chla.sym),
            sym.cm2= mean(sym.cm2),
            prot_mg.cm2= mean( prot_mg.cm2),
            Host_AFDW_mg.cm2= mean(Host_AFDW_mg.cm2),
            Sym_AFDW_mg.cm2= mean(Sym_AFDW_mg.cm2),
            CM2.year= mean(CM2.year))
            
# make new column that will have final treatment 
raw$final <- ifelse(
  raw$treatment %in% c("PS", "SS"), "5", 
  ifelse(raw$treatment %in% c("PP", "SP"), "18", NA))
raw$home_away <- ifelse(
  raw$treatment %in% c("PP", "SS"), "Home", 
  ifelse(raw$treatment %in% c("PS", "SP"), "Away", NA))

# calculate adult calices 
raw <- raw %>%
  mutate(Dadult = (D2 - Dsmall))

# make new column that will have final treatment 
PAM$final <- ifelse(
  PAM$full_treatment %in% c("OFAV_PS", "OFAV_SS", "OFRA_PS"), "5", 
  ifelse(PAM$full_treatment %in% c("OFAV_PP", "OFAV_SP", "OFRA_PP"), "18", NA))
PAM$home_away <- ifelse(
  PAM$full_treatment %in% c("OFAV_PP", "OFRA_PP", "OFAV_SS"), "Home", 
  ifelse(PAM$full_treatment %in% c("OFRA_PS", "OFAV_SP", "OFAV_PS"), "Away", NA))

# Define my groups for stat compare means 
treatment_comparisonsx <- list( c("OFAV_PS","OFAV_PP"), c("OFAV_SP","OFAV_SS"), c("OFRA_PP","OFRA_PS"), c("OFRA_PP","OFAV_PP"), c("OFAV_PP","OFAV_SS"))
treatment_comparisons <- list( c("OFAV_PS","OFAV_PP"), c("OFAV_SP","OFAV_SS"), c("OFRA_PP","OFRA_PS"))

# set orders for graphs 
measurement_order <- c('OFAV_S','OFAV_P','OFRA_P') 
x_order <- c('OFAV_SS','OFAV_SP','OFAV_PS','OFAV_PP','OFRA_PS','OFRA_PP') 

# define offset of points 
pj <- position_jitterdodge(jitter.width=0.4, seed=9, jitter.height = 0)

# NOTES: 
# Remove outliers: Each variable in this data set has been filtered for outliers > 3x z-score 
# Check for Normalcy: Run Shapiro and Bartlett tests on all parameters to test for normality and homogeneity 
# All physiology metrics do not have parametric shapes or homogeneity, must use non-parametric tests on all analyses 
# Use Wilcoxon's tests for everything 

blank <- ggplot(raw, aes(x=final, y = D2, fill = home_away)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#FF6347", "#A3A3D9"), labels = c("Away", "Home")) + 
  labs(y= "", x = "", fill='Treatment') + 
  theme_classic() +
  theme(legend.position=c(0.35, 0.3), 
        legend.background = element_rect(fill = "white", color = NA), # White background
        legend.margin = margin(200, 150, 400, 150),
        text = element_text(size=40), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        legend.key.size = unit(3, "cm"),       # Increase overall size of the legend keys
        plot.margin=unit(c(1,0,0,1),"cm")
        )

blank
ggsave("TL_TRANS_blank.jpg", plot = blank, path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 5, height = 5)

# create custom color pallette 
custom_colors <- c("Home" = "#A3A3D9", "Away" = "#FF6347")  # Adjust as needed

# Reaction Norm Calculation -----------------------------------------------

# combine RAW and PAM values 
PAM <- PAM %>% mutate(colony_id = as.character(colony_id))
full_phys <- full_join(raw, PAM)
full_phys <- full_phys %>%
  select(species, group, full_treatment, final, home_away, CM2.year, D2, Dsmall, Dadult, A, CA, di, Cdi, Qm, Fv.Fm, prot_mg.cm2, sym.cm2, Sym_AFDW_mg.cm2, Host_AFDW_mg.cm2, chla.ug.cm2)

phys_longer <- full_phys %>%
  pivot_longer(cols = c("A", "CA", "di", "Cdi", "D2", "Dsmall", "Dadult", "chla.ug.cm2", "sym.cm2", "prot_mg.cm2", "Host_AFDW_mg.cm2", "Sym_AFDW_mg.cm2", "CM2.year", "Qm", "Fv.Fm"), names_to = "metric", values_to = "value") %>%
  filter(value != "NA") %>%
  select(group, final, home_away, metric, value)

reaction_norms <- phys_longer %>%
  group_by(group, final, metric, home_away)%>%
  summarise(mean = mean(value),
            sd = sd(value))

# Fig 2. Environmental Data ------------------------------------------------------

# prep temp & DO data 
temp_DO <- temp_DO %>% 
  filter(!is.na(datetime)) %>% 
  mutate(datetime = mdy_hm(datetime)) 
temp_DO$date <- as.Date(temp_DO$datetime, "%Y-%m-%d %H:%M:%S")

# prep light data
light_shallow$date <- as.Date(light_shallow$date, "%d/%m/%y")
light_deep$date <- as.Date(light_deep$date, "%d/%m/%y")
  
light <- rbind(light_shallow, light_deep)
light$datetime <- as.POSIXct(paste(light$date, light$time), format = "%Y-%m-%d %H:%M:%S")

# Check data 
ggplot(light, aes(x=datetime, y=LightRaw, color = depth)) + geom_point()

# Combine Files 
enviro_merged <- merge(light, temp_DO, all=TRUE) %>%
  select(datetime, date, time, depth, LightRaw, DO_mg.L, Temp_C)

write_csv(enviro_merged, '~/Desktop/GitHub/TLPR21/Meta_Environmental_Data/TL_TRANS_Enviro_all_UPDATE.csv')

# include only times 1hr before and after sunrise (6am) and set (7pm)

# filter by time of day 
enviro_merged$time <- as_hms(enviro_merged$time)
enviro_merged$daylight <- ifelse(enviro_merged$time >= hms::as_hms('5:00:00') & enviro_merged$time <= hms::as_hms('19:00:00'), paste(enviro_merged$LightRaw), NA) %>%
  as.numeric(enviro_merged$daylight)  

# plot by time of day 
ggplot(enviro_merged, aes(x=time, y=daylight)) + geom_point()

# Daily Means Graphs

#merged <- read.csv('~/Desktop/GITHUB/TLPR21/Meta_Environmental_Data/TL_TRANS_Enviro_all.csv') 
#merged$datetime <- as.POSIXct(merged$datetime, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")

# daily means 
daily_light <- enviro_merged %>%
  group_by(date,depth) %>%
  summarize(mean_light = mean(LightRaw, na.rm = TRUE),
            sd_light = sd(LightRaw, na.rm = TRUE)) %>%
  ungroup(date) %>%
  filter(.$date < as.Date("2021-08-14")) %>%
  mutate(sd_min = mean_light-sd_light) %>%
  mutate(sd_max = mean_light+sd_light)

daily_DO <- enviro_merged %>%
  group_by(date,depth) %>%
  summarize(mean_DO = mean(DO_mg.L, na.rm = TRUE), 
            sd_DO = sd(DO_mg.L, na.rm = TRUE)) %>%
  mutate(sd_min = mean_DO-sd_DO) %>%
  mutate(sd_max = mean_DO+sd_DO)

daily_temp <- enviro_merged %>%
  group_by(date,depth) %>%
  summarize(mean_temp = mean(Temp_C, na.rm = TRUE),
            sd_temp = sd(Temp_C, na.rm = TRUE))%>%
  mutate(sd_min = mean_temp-sd_temp) %>%
  mutate(sd_max = mean_temp+sd_temp)

# GRAPHING DAILY VALUES 

# TEMP
enviro_temp <- ggplot(daily_temp, aes(x=date, y=mean_temp, color=depth)) + 
  geom_ribbon(aes(ymin=sd_min, ymax=sd_max, color=depth, fill = depth), alpha = 0.3, colour = NA) + ###
  geom_line(size=2) +
  geom_point(size=4) + 
  scale_color_manual(values = c("#3E3ED1", "#66BBBB")) + 
  scale_fill_manual(values = c("#3E3ED1", "#66BBBB")) +
  # Aesthetics 
  theme_bw() + 
  theme(
    legend.position="none", 
    text = element_text( size=35), 
    axis.text.x = element_blank(), 
    axis.ticks.x= element_blank()
  )+ 
  labs(y = "Mean Temperature (˚C)", x = "") +
  scale_x_date(limits = as.Date(c('2021-07-12','2021-09-01')))
enviro_temp

# Light
enviro_light <- ggplot(daily_light, aes(x=date, y=mean_light, color=depth)) + 
  geom_ribbon(aes(ymin=sd_min, ymax=sd_max, color=depth, fill = depth), alpha = 0.3, colour = NA) + ###
  geom_line(size=2) + 
  geom_point(size=4) + 
  scale_color_manual(values = c("#3E3ED1", "#66BBBB"), labels=c('18m', '5m'), name="Depth") +
  scale_fill_manual(values = c("#3E3ED1", "#66BBBB"), labels=c('18m', '5m'), name="Depth") +
  # Aesthetics 
  theme_bw() + 
  theme(
    legend.position=c(.8, .5), 
    text = element_text( size=35), 
    legend.box.background = element_rect(color="black", size=1.5), 
    axis.text.x = element_blank(), 
    axis.ticks.x= element_blank()
  ) +
  labs(y = expression("Mean Daily Light (lum ft"^{-2}*")"), x = "") +
  scale_x_date(limits = as.Date(c('2021-07-12','2021-09-01')))
enviro_light

# DO 
enviro_DO <- ggplot(daily_DO, aes(x=date, y=mean_DO, color=depth)) + 
  geom_ribbon(aes(ymin=sd_min, ymax=sd_max, color=depth, fill = depth), alpha = 0.3, colour = NA) +
  geom_line(size=2) + 
  geom_point(size=4) + 
  scale_color_manual(values = c("#3E3ED1", "#66BBBB")) +
  scale_fill_manual(values = c("#3E3ED1", "#66BBBB")) +
  # Aesthetics 
  theme_bw() +
  theme(
    legend.position="none", 
    text = element_text(size=35), 
  ) + 
  labs(y = "Mean Dissolved Oxygen (mg/L)", x = "Date") +
  scale_x_date(limits = as.Date(c('2021-07-12','2021-09-01')))

#combine ecotype and treatment 
enviro_arrange <- plot_grid(enviro_temp, enviro_light, enviro_DO, 
                      ncol = 1, align = "v",
                      labels = c("A", "B", "C"),  label_size = 35, label_x = 0.11, label_y = 0.99)
enviro_arrange

ggsave("TL_TRANS_Enviro_Arrange.jpg", plot = enviro_arrange, path = '~/Desktop/TL_TRANS_Manuscript/Graphs/', width = 15, height = 25)

# Compare shallow and deep means: 
merged_longer <- enviro_merged %>%
  pivot_longer(!c(datetime, date, time, depth), names_to = "measurement", values_to = "value")
ggplot(merged_longer, aes(x=depth, y=value)) + 
  geom_boxplot() + 
  stat_compare_means(method = "t.test") +
  facet_wrap(~measurement,scales = "free")

# Compare means & write df 
enviro_means <- merged_longer %>%
  group_by(depth, measurement) %>%
  summarise(mean = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE))
write.csv(enviro_means , "~/Desktop/TL_TRANS_Manuscript/STATS/TL_TRANS_STATS_Enviro_means.csv", row.names = FALSE)

# Paired t-test 
# pivot wider 
enviro_wider <- enviro_merged %>%
  group_by(datetime) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from=depth, values_from = c(LightRaw, DO_mg.L, Temp_C))

# function 
perform_t_test <- function(group1, group2, group_label) {
  t_test <- t.test(group1, group2, paired = FALSE)
  data.frame(
    Comparison = group_label,
    t = t_test$statistic,
    df = t_test$parameter,
    p = formatC(t_test$p.value, format = "e", digits = 3)
  )
}

# Perform the t-tests and collect results
test_results_list <- list(
  perform_t_test(enviro_wider$Temp_C_shallow, enviro_wider$Temp_C_deep, "temp"),
  perform_t_test(enviro_wider$LightRaw_shallow, enviro_wider$LightRaw_deep, "Light"),
  perform_t_test(enviro_wider$DO_mg.L_shallow, enviro_wider$DO_mg.L_deep, "DO")
)

# Combine all results into one data frame
all_enviro <- do.call(rbind, test_results_list)

write.csv(all_enviro, "~/Desktop/TL_TRANS_Manuscript/STATS/TL_TRANS_STATS_Enviro_t-test.csv", row.names = FALSE)

# **Summary of environmental data**
# Temperature: shallow is significantly higher by 0.02 deg C 
# DO: shallow is significantly lower by 0.44 mg/L DO
# Light (all): shallow is significantly higher by 118.44
# Light (Daylight only): Shallow is significantly higher by 70 lum

# Fig 3a. Survival Curves  ----------------------------------------------------------------

# filter out samples that went missing 
surv_filtered <- surv %>% filter(!is.na(date.planted)) %>% filter(!is.na(date.dead)) %>% filter(!is.na(final.status))

#reformat all dates to POSIXt format using lubridate and delete old columns
surv_filtered$ymd.planted <- ymd(surv_filtered$date.planted)
surv_filtered$ymd.mortality <- ymd(surv_filtered$date.dead)
surv_new <- select(surv_filtered, -date.planted, -date.dead)

# earliest and latest days 
surv_check <- surv_new %>%
  group_by(full_treatment) %>%
  filter(!is.na(ymd.planted)) %>%
  filter(!is.na(ymd.mortality)) %>%
  summarise(plant_early = min(ymd.planted),
            plant_late = max(ymd.planted),
            end_early = min(ymd.mortality),
            end_late = max(ymd.mortality),
            max_days = as.numeric(plant_early - end_late))

#calculate number of days between planting and mortality
surv_new$days.alive <- difftime(surv_new$ymd.mortality , surv_new$ymd.planted, units = c("days"))

# calculate n alive and dead and % mortality 
n_survival <- surv_filtered %>%
  group_by(final.status, group, final) %>%
  summarise(count = n())  %>% 
  pivot_wider(names_from = final.status, values_from = count) %>%
  replace(is.na(.), 0)  %>%
  mutate(total = `0` + `1`) %>%
  mutate(mean = (`0`/total)*100) #this is percent survival

#make variables numeric 
surv_new$days.alive <- as.numeric(surv_new$days.alive)
surv_new$final.status <- as.numeric(surv_new$final.status)

# create a surv object 
all_treatments <- Surv(surv_new$days.alive, surv_new$final.status)
fit1 <- survfit(all_treatments ~ full_treatment, data = surv_new)
plot_treatment <- ggsurvplot(fit1, data = surv_new, pval = TRUE, legend = "bottom", legend.title = "Treatment")
plot_treatment <- plot_treatment + xlab("Days")
plot_treatment

### OFAV S 

# plot by treatment- comparisons 
new_OFAVS <- surv_new %>% filter(full_treatment == "OFAV_SP" | full_treatment == "OFAV_SS")

### Remove 0 days alive 
new_OFAVS <- new_OFAVS %>% filter(days.alive != 0)

# fit the K-M model 
OFAVS <- Surv(new_OFAVS$days.alive, new_OFAVS$final.status)
fit_OFAVS <- survfit(OFAVS ~ full_treatment, data = new_OFAVS)

# PLOT OFAV Shallow COMPARISON 
plot_OFAVS  <- ggsurvplot(fit_OFAVS , data = new_OFAVS, 
                          pval = FALSE, 
                          #pval.coord = c(28, 0.24),
                          size=1.5,
                          palette = c("#FF6347", "#A3A3D9"),
                          ylim = c(0.1, 1),
                          xlim = c(0, 128),
                          ylab= c("Survival Probability"), 
                          xlab = c(""),
                          ggtheme = theme_classic2(base_size=30),
                          conf.int = 0.95,
                          break.time.by = 20)

plot_OFAVS$plot <- plot_OFAVS$plot +  
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 10), 
        legend.position = "none",  # Remove legend
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"))   # Bold legend title)

plot_OFAVS

### OFAV DEEP

# plot by treatment- comparisons 
new_OFAVP <- surv_new %>% filter(full_treatment == "OFAV_PP" | full_treatment == "OFAV_PS")

### Remove 0 days alive 
new_OFAVP <- new_OFAVP %>% filter(days.alive != 0)

# Fit the K-M model 
OFAVP <- Surv(new_OFAVP$days.alive, new_OFAVP$final.status)
fit_OFAVP <- survfit(OFAVP ~ full_treatment, data = new_OFAVP)

# PLOT OFAV DEEP COMPARISON 
plot_OFAVP  <- ggsurvplot(fit_OFAVP , data = new_OFAVP, 
                          pval = FALSE, 
                          #pval.coord = c(.3, 0.8),
                          size=1.5,
                          palette = c("#A3A3D9", "#FF6347"),
                          ylim = c(0.1, 1),
                          #xlim = c(-105, 23),
                          xlim = c(0, 23),
                          xlab= c(""),
                          ylab= c(""), 
                          ggtheme = theme_classic2(base_size=30),
                          conf.int = 0.95,
                          break.time.by = 10)

plot_OFAVP$plot <- plot_OFAVP$plot + 
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 10), 
        legend.position = "none",  # Remove legend
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),   # Remove y-axis label
        axis.text.y = element_blank(),    # Remove y-axis numbers
        axis.ticks.y = element_blank())   # Bold legend title
plot_OFAVP


### OFRA DEEP

# plot by treatment- comparisons 
new_OFRA <- surv_new %>% filter(full_treatment == "OFRA_PP" | full_treatment == "OFRA_PS")

### Remove 0 days alive 
new_OFRA <- new_OFRA %>% filter(days.alive != 0)

# fit the K-M model 
OFRA <- Surv(new_OFRA$days.alive, new_OFRA$final.status)
fit_OFRA <- survfit(OFRA ~ full_treatment, data = new_OFRA)
# PLOT OFRA COMPARISON 
plot_OFRA  <- ggsurvplot(fit_OFRA, data = new_OFRA, 
                             pval = FALSE, 
                             #pval.coord = c(.2, 0.95),
                             size=1.5,
                             legend.title = "Treatment Category", 
                             legend.labs = c("Home", "Away"), 
                             palette = c("#A3A3D9", "#FF6347"), 
                             ylim = c(0.1, 1),
                             xlim = c(0, 73),
                            #xlim = c(-55, 73),
                             ylab = c(""),
                             xlab= c(""),
                             ggtheme = theme_classic2(base_size=30),
                             conf.int = 0.95, 
                             break.time.by = 20)

plot_OFRA$plot <- plot_OFRA$plot +  
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 10), 
        scale_x_continuous(breaks = seq(-55, 73, by = 10)),   # Set tick marks every 25 units
        #axis.title.x = element_text(face = "bold", size = 20),
       # legend.title =  element_text(face = "bold", size = 30),
       # legend.position = "right",
       legend.position = "none",  # Remove legend
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),   # Remove y-axis label
        axis.text.y = element_blank(),    # Remove y-axis numbers
        axis.ticks.y = element_blank())
plot_OFRA

# P-values: YELLOW: p<0.0001, PINK: 0.0049, PURPLE: 0.00017

# Convert ggsurvplot objects to grobs
grob_OFAVS <- ggplotGrob(plot_OFAVS$plot)
grob_OFAVP <- ggplotGrob(plot_OFAVP$plot)
grob_OFRA <- ggplotGrob(plot_OFRA$plot)

# Combine the plots into a single layout
survival_arrange1 <- grid.arrange(
  grob_OFAVS, grob_OFAVP, grob_OFRA,  # Use the grob versions
  nrow = 1,
  widths = c(0.63, 0.10, 0.27),
  bottom = text_grob("Days Since Transplantation", size = 30, face = "bold")
  )  

#save graphs 
ggsave("TL_TRANS_Surival_Arrange_new.jpg", plot = survival_arrange1 , path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 18, height = 5)

# SURVIVAL STATS 

# Define a function to perform a log-rank test and extract the results
perform_log_rank_test <- function(surv_object, grouping_variable, data, x) {
  # Perform the log-rank test
  log_rank_test <- survdiff(surv_object ~ get(grouping_variable), data = data)
  
  # Extract test statistic, degrees of freedom, and p-value
  test_statistic <- log_rank_test$chisq
  degrees_of_freedom <- length(log_rank_test$n) - 1
  p_value <- 1 - pchisq(test_statistic, df = degrees_of_freedom)
  
  # Return a data frame with the results
  data.frame(
    comparison = x,
    chisq = test_statistic,
    df = degrees_of_freedom,
    p = formatC(p_value, digits = 3)
  )
}

# Example: List of survival objects and grouping variables
tests <- list(
  list(surv_object = OFAVP, grouping_variable = "full_treatment", data = new_OFAVP, x = "OFAV Deep"),
  list(surv_object = OFAVS, grouping_variable = "full_treatment", data = new_OFAVS, x = "OFAV Shallow"),
  list(surv_object = OFRA, grouping_variable = "full_treatment", data = new_OFRA, x = "OFRA")
)

# Apply the function to all tests
all_log_rank_results <- do.call(rbind, lapply(tests, function(test) {
  perform_log_rank_test(test$surv_object, test$grouping_variable, test$data, test$x)
}))

# Save 
write.csv(all_log_rank_results, "~/Desktop/TL_TRANS_Manuscript/STATS/TL_TRANS_STATS_survival_kaplan-meier.csv", row.names = FALSE)

# FIG Sx. Parametric survival curves  ---------------------------------------------

# using the exponential shape, we assume that rate of survivorship is fixed over time 

#fit1_OFAVS
weib_OFAVS <- survreg(Surv(days.alive, final.status) ~ full_treatment, data = new_OFAVS, dist = "weibull")
summary(weib_OFAVS) # treatment p < 2e-16
extractAIC(weib_OFAVS) # Weibull df = 3, AIC =2667.557

exp_OFAVS <- survreg(Surv(days.alive, final.status) ~ full_treatment, data = new_OFAVS, dist = "exp")
summary(exp_OFAVS) # treatment p <2e-16
extractAIC(exp_OFAVS) # exp df = 2, AIC = 2728.209

#fit1_OFAVP
weib_OFAVP <- survreg(Surv(days.alive, final.status) ~ full_treatment, data = new_OFAVP, dist = "weibull")
summary(weib_OFAVP) # treatment p = 0.014
extractAIC(weib_OFAVP) # Weibull df = 3, AIC = 1380.831

exp_OFAVP <- survreg(Surv(days.alive, final.status) ~ full_treatment, data = new_OFAVP, dist = "exp")
summary(exp_OFAVP) # treatment p = 0.0078
extractAIC(exp_OFAVP) # exp df = 2, AIC = 1528.609

#fit1_OFRA
weib_OFRA <- survreg(Surv(days.alive, final.status) ~ full_treatment, data = new_OFRA, dist = "weibull")
summary(weib_OFRA) # treatment p = 0.997
extractAIC(weib_OFRA) # Weibull df = 3, AIC = 275.2446

exp_OFRA <- survreg(Surv(days.alive, final.status) ~ full_treatment, data = new_OFRA, dist = "exp")
summary(exp_OFRA) # treatment p = 1
extractAIC(exp_OFRA) # exp df = 2, AIC = 281.79

# In all cases, weibull has a lower AIC than exponential. 
# This allows the equation to assume that the instantaneous likelyhood of survival changes over the course of the acclimation period 

# OFAVS
weib.treat.OFAVSS = predict(weib_OFAVS, newdata=list(full_treatment="OFAV_SS"),type="quantile",p=seq(.01,.99,by=.01))
weib.treat.OFAVSP = predict(weib_OFAVS, newdata=list(full_treatment="OFAV_SP"),type="quantile",p=seq(.01,.99,by=.01))
exp.treat.OFAVSS = predict(exp_OFAVS, newdata=list(full_treatment="OFAV_SS"),type="quantile",p=seq(.01,.99,by=.01))
exp.treat.OFAVSP = predict(exp_OFAVS, newdata=list(full_treatment="OFAV_SP"),type="quantile",p=seq(.01,.99,by=.01))

# OFAVP
weib.treat.OFAVPP = predict(weib_OFAVP, newdata=list(full_treatment="OFAV_PP"),type="quantile",p=seq(.01,.99,by=.01))
weib.treat.OFAVPS = predict(weib_OFAVP, newdata=list(full_treatment="OFAV_PS"),type="quantile",p=seq(.01,.99,by=.01))
exp.treat.OFAVPP = predict(exp_OFAVP, newdata=list(full_treatment="OFAV_PP"),type="quantile",p=seq(.01,.99,by=.01))
exp.treat.OFAVPS = predict(exp_OFAVP, newdata=list(full_treatment="OFAV_PS"),type="quantile",p=seq(.01,.99,by=.01))

# OFRA
weib.treat.OFRAPP = predict(weib_OFRA, newdata=list(full_treatment="OFRA_PP"),type="quantile",p=seq(.01,.99,by=.01))
weib.treat.OFRAPS = predict(weib_OFRA, newdata=list(full_treatment="OFRA_PS"),type="quantile",p=seq(.01,.99,by=.01))
exp.treat.OFRAPP = predict(exp_OFRA, newdata=list(full_treatment="OFRA_PP"),type="quantile",p=seq(.01,.99,by=.01))
exp.treat.OFRAPS = predict(exp_OFRA, newdata=list(full_treatment="OFRA_PS"),type="quantile",p=seq(.01,.99,by=.01))

#### Add the curves to the plots 

#OFAVS
df_OFAVS_w = data.frame(y=seq(.99,.01,by=-.01), treat1=weib.treat.OFAVSS, treat2=weib.treat.OFAVSP)
OFAVS_long_weib = gather(df_OFAVS_w, key= "treat", value="time", -y)
df_OFAVS_e = data.frame(y=seq(.99,.01,by=-.01), treat1=exp.treat.OFAVSS, treat2=exp.treat.OFAVSP)
OFAVS_long_exp = gather(df_OFAVS_e, key= "treat", value="time", -y)

plot_OFAVS$plot = plot_OFAVS$plot + 
  geom_line(data=OFAVS_long_weib, aes(x=time, y=y, group=treat), color = "gold", size = 0.8) +
  geom_line(data=OFAVS_long_exp, aes(x=time, y=y, group=treat), color="darkblue", size = 0.8) 
plot_OFAVS

#OFAVS
df_OFAVP_w = data.frame(y=seq(.99,.01,by=-.01), treat1=weib.treat.OFAVPP, treat2=weib.treat.OFAVPS)
OFAVP_long_weib = gather(df_OFAVP_w, key= "treat", value="time", -y)
df_OFAVP_e = data.frame(y=seq(.99,.01,by=-.01), treat1=exp.treat.OFAVPP, treat2=exp.treat.OFAVPS)
OFAVP_long_exp = gather(df_OFAVP_e, key= "treat", value="time", -y)

plot_OFAVP$plot = plot_OFAVP$plot + 
  geom_line(data=OFAVP_long_weib, aes(x=time, y=y, group=treat), color = "gold", size = 0.8) +
  geom_line(data=OFAVP_long_exp, aes(x=time, y=y, group=treat), color="darkblue", size = 0.8) 
plot_OFAVP

#OFAVS
df_OFRA_w = data.frame(y=seq(.99,.01,by=-.01), treat1=weib.treat.OFRAPP, treat2=weib.treat.OFRAPS)
OFRA_long_weib = gather(df_OFRA_w, key= "treat", value="time", -y)
df_OFRA_e = data.frame(y=seq(.99,.01,by=-.01), treat1=exp.treat.OFRAPP, treat2=exp.treat.OFRAPS)
OFRA_long_exp = gather(df_OFRA_e, key= "treat", value="time", -y)

plot_OFRA$plot = plot_OFRA$plot + 
  geom_line(data=OFRA_long_weib, aes(x=time, y=y, group=treat), color = "gold", size = 0.8) +
  geom_line(data=OFRA_long_exp, aes(x=time, y=y, group=treat), color="darkblue", size = 0.8) 
plot_OFRA

#### 
grob_OFAVS <- ggplotGrob(plot_OFAVS$plot)
grob_OFAVP <- ggplotGrob(plot_OFAVP$plot)
grob_OFRA <- ggplotGrob(plot_OFRA$plot)

# Combine the plots into a single layout
para_survival_arrange <- grid.arrange(
  arrangeGrob(grob_OFAVS, top = textGrob("OFAV Shallow", gp = gpar(fontsize = 25, fontface="bold")), padding = unit(1, "lines")),
  arrangeGrob(grob_OFAVP, top = textGrob("OFAV Deep", gp = gpar(fontsize = 25, fontface="bold")), padding = unit(1, "lines")),
  arrangeGrob(grob_OFRA, top = textGrob("OFRA Deep", gp = gpar(fontsize = 25, fontface="bold")), padding = unit(1, "lines")),
  nrow = 1,
  widths = c(0.63, 0.10, 0.27),
  bottom = text_grob("Days Since Transplantation", size = 30, face = "bold")
)  

#save graphs
ggsave("TL_TRANS_Surival_Arrange_parametric_new.jpg", plot = para_survival_arrange , path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 15, height = 6)

#

# Fig 3b. Growth rate -------------------------------------------------------------

# pull mean growth rates 
mean_growth <- reaction_norms %>%
  filter(metric == "CM2.year")

### GROWTH BOXPLOT 

# factor 
raw$group <- factor(raw$group, levels = c("OFAV_S", "OFAV_P", "OFRA_P")) 
raw$final <- factor(raw$final, levels = c("5", "18")) 

# growth boxplots 
box_growth <- ggplot(raw, aes(x=final, y=CM2.year)) +
  geom_boxplot(aes(fill=home_away), alpha=0.7, outlier.size=0) +
  geom_point(aes(color=home_away), position = pj, size=3) + 
  facet_wrap(~ group, strip.position = "top", labeller = labeller(group = custom_group_labels)) +  # Move facet labels below
  theme_bw() +
  # LABELS & AXES 
  labs(x = "Transplant Depth", y = "Growth (cm²/year)", color = "Habitat") +
  scale_color_manual(values = c("black", "black")) +  # Set custom point colors
  scale_fill_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2() + 
  theme(text = element_text(size=30), 
        #strip.text = element_text(size = 14, face = "bold"),
        #legend.position = c(0.85, 0.9),  # Move legend inside (adjust x,y for fine-tuning)
        #legend.background = element_rect(fill = alpha("white", 0.7)),  # Transparent background
        #legend.key = element_blank(),
        legend.position = "none",  # Remove legend
        axis.title.x =  element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin=unit(c(0,0,0,2),"cm"),
        axis.line.y = element_line(colour = "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5)
  
ggsave("TL_TRANS_growth_box.jpg", plot = box_growth, path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 20, height = 15)

### GROWTH RACTION NORM

# factor 
mean_growth$group <- factor(mean_growth$group, levels = c("OFAV_S", "OFAV_P", "OFRA_P")) 
mean_growth$final <- factor(mean_growth$final, levels = c("5", "18")) 
mean_growth$home_away<- factor(mean_growth$home_away, levels = c("Home", "Away")) 

reaction_norm_growth <- ggplot(mean_growth, aes(x=final, y=mean)) +
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_blank(), 
        legend.position = "none",  # Remove legend
        axis.title.x =  element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0,0,0,2),"cm")
  ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) + 
 # scale_y_continuous(labels = scales::number_format(accuracy = 0.)) + 
  labs(x = "Depth (m)", y = "Growth (cm²/year)")
reaction_norm_growth

# Fig 3d. Calice types ------------------------------------------------------------

buds_means <- raw %>%
  group_by(full_treatment, group, final, home_away) %>%
  summarise(
    mean_bud = mean(Dsmall, na.rm = TRUE), 
    mean_adult = mean(Dadult, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(mean_adult, mean_bud), names_to = "type", values_to = "Value")

# Ensure mean_bud appears above mean_adult
buds_means$type <- factor(buds_means$type, levels = c("mean_bud", "mean_adult")) 
buds_means$group <- factor(buds_means$group, levels = c("OFAV_S", "OFAV_P", "OFRA_P")) 
buds_means$final <- factor(buds_means$final, levels = c("5", "18")) 

# Labels
custom_group_labels <- c("OFAV_S" = "OFAV Shallow", "OFAV_P" = "OFAV Deep", "OFRA_P" = "OFRA Deep")

bar_buds <- ggplot(buds_means, aes(x = final, y = Value, fill = home_away, alpha = type)) +
  # DATA with pattern
  geom_bar_pattern(stat = "identity", 
                   aes(pattern = type),  # Map pattern to type
                   pattern_density = 0.3,  # Adjust density of polka dots
                   pattern_spacing = 0.05,  # Adjust spacing
                   pattern_fill = "white",  # Color of polka dots
                   pattern_color = "white") +  # Outline color
  # AESTHETICS 
  scale_fill_manual(values = custom_colors,  
                    labels = c("mean_adult" = "Adults", "mean_bud" = "Buds")) +  
  scale_alpha_manual(values = c("mean_adult" = 1, "mean_bud" = 0.6)) +  
  scale_pattern_manual(values = c("mean_adult" = "none", "mean_bud" = "circle")) +  # Define patterns
  facet_wrap(~ group, strip.position = "top", labeller = labeller(group = custom_group_labels)) +  
  scale_y_continuous(limits = c(0,4), expand = expansion(mult = c(0, 0.05))) + 
  # LABELS & AXES 
  labs(x = "Transplant Depth (m)", y = "Calice density (per cm²)", fill = "Calice Type") +
  theme_classic2() + 
  theme(text = element_text(size=30), 
        legend.position = "none", 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        #legend.position = c(0.85, 0.9),  
        #legend.background = element_rect(fill = alpha("white", 0.7)),  
        #legend.key = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin=unit(c(0,0,0,2),"cm"))
       

bar_buds

ggsave("TL_TRANS_buds.jpg", plot = bar_buds, path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 25, height = 20)

# Fig 3. Fitness Arrange -----------------------------------------------

### STACKED ARANGEMENT 

fitness <- grid.arrange(
  arrangeGrob(survival_arrange1, top = textGrob("", gp = gpar(fontsize = 10)), padding = unit(3, "lines")),
  arrangeGrob(reaction_norm_growth, top = textGrob("", gp = gpar(fontsize = 10)), padding = unit(6, "lines")), 
  arrangeGrob(bar_buds, top = textGrob("", gp = gpar(fontsize = 10)), padding = unit(2, "lines")),
  nrow = 3,
  bottom = grobTree(
    textGrob("A                               OFAV Shallow                     OFAV Deep    OFRA Deep", x = 0.01, y = 245, just = "left", gp = gpar(fontsize = 40, fontface = "bold")),
    textGrob("B            OFAV Shallow                    OFAV Deep                    OFRA Deep", x = 0.01, y = 160, just = "left", gp = gpar(fontsize = 40, fontface = "bold")),
    textGrob("C", x = 0.01, y = 80, just = "left", gp = gpar(fontsize = 40, fontface = "bold")),
    textGrob("*", x = 0.38, y = 210, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("*", x = 0.68, y = 210, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("*", x = 0.86, y = 210, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("*", x = 0.22, y = 135, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("*", x = 0.19, y = 71, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("bud", x = 0.22, y = 73, just = "left", gp = gpar(fontsize = 30)),
    textGrob("*", x = 0.46, y = 71, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("bud, adult, all", x = 0.49, y = 73, just = "left", gp = gpar(fontsize = 30)),
    textGrob("*", x = 0.8, y = 71, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("bud, all", x = 0.83, y = 73, just = "left", gp = gpar(fontsize = 30))
    ))

#ggsave("TL_TRANS_fitness.jpg", plot = fitness, path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 20, height = 25)

### GRID ARANGEMENT 

fit1 <- grid.arrange( 
  arrangeGrob(survival_arrange1, top = textGrob("", gp = gpar(fontsize = 10)), padding = unit(3, "lines")), 
  arrangeGrob(survival_arrange1, top = textGrob("", gp = gpar(fontsize = 10)), padding = unit(3, "lines")), 
  nrow = 2, heights = c(0.56, 0.44) )

fit2 <- grid.arrange(
  arrangeGrob(reaction_norm_growth, top = textGrob("", gp = gpar(fontsize = 10)), padding = unit(3, "lines")), 
  arrangeGrob(bar_buds, top = textGrob("", gp = gpar(fontsize = 10)), padding = unit(3.5, "lines")),
  nrow = 2, heights = c(0.5, 0.5) )

fitness <- grid.arrange(
 fit1, fit2,
  nrow = 1, 
  bottom = grobTree(
    textGrob("                           OFAV Shallow           OFAV Deep   OFRA Deep                OFAV Shallow           OFAV Deep            OFRA Deep", x = 0.01, y = 145, just = "left", gp = gpar(fontsize = 36, fontface = "bold")),
    textGrob("A", x = 0.01, y = 145, just = "left", gp = gpar(fontsize = 40, fontface = "bold")),
    textGrob("B", x = 0.52, y = 145, just = "left", gp = gpar(fontsize = 40, fontface = "bold")),
    textGrob("D", x = 0.52, y = 70, just = "left", gp = gpar(fontsize = 40, fontface = "bold")),
    
    textGrob("*", x = 0.19, y = 112, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("*", x = 0.34, y = 112, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("*", x = 0.43, y = 112, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("*", x = 0.63, y = 112, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    
    textGrob("*", x = 0.6, y = 63, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("bud", x = 0.62, y = 65, just = "left", gp = gpar(fontsize = 28)),
    textGrob("*", x = 0.72, y = 63, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("bud, adult, all", x = 0.74, y = 65, just = "left", gp = gpar(fontsize = 28)),
    textGrob("*", x = 0.9, y = 63, just = "left", gp = gpar(fontsize = 80, fontface = "bold")),
    textGrob("bud, all", x = 0.92, y = 65, just = "left", gp = gpar(fontsize = 28))
  ))

ggsave("TL_TRANS_fitness.jpg", plot = fitness, path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 30, height = 15)

# Fig 4. Morph Reaction Norms --------------------------------------------------------

# Assign orders for each variable 
reaction_norms$group <- factor(reaction_norms$group, levels = c("OFAV_S", "OFAV_P", "OFRA_P")) 
reaction_norms$home_away <- factor(reaction_norms$home_away, levels = c("Home", "Away"))
reaction_norms$final <- factor(reaction_norms$final, levels = c("5", "18")) 

facet_labels <- c(
  "OFAV_S" = "OFAV Shallow", 
  "OFAV_P" = "OFAV Deep",
  "OFRA_P" = "OFRA Deep"
)

morph_norms <- reaction_norms %>%
  filter(metric %in% c("A"))

reaction_norm_plot <- ggplot(morph_norms, aes(x=final, y=mean)) +
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(strip.placement = "outside", 
        text = element_text(size=30),
        legend.position=c(0.3, 0.8), 
        legend.background = element_rect(fill = "white", color = NA),
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_text(face = "bold"), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        ) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(y = "Corallite Area (cm²)", color = "Treatment")
reaction_norm_plot

morph_norms2 <- reaction_norms %>%
  filter(metric %in% c("CA"))

reaction_norm_plot2 <- ggplot(morph_norms2, aes(x=final, y=mean)) +
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        legend.position = "none", 
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_blank(), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold")
  ) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Depth (m)", y = "Calice Area (cm²)")
reaction_norm_plot2


morph2 <- grid.arrange(reaction_norm_plot, reaction_norm_plot2, nrow=2,
                       bottom = grobTree(
                         textGrob("*", x = 0.83, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold"))
                       ))

ggsave("TL_TRANS_reaction_norms_morph.jpg", plot = morph2, path = '~/Desktop/TL_TRANS_Manuscript/Graphs/', width = 9, height = 18)

# Fig 5. Reaction Norms phys -----------------------------------------------------

### SYM 
norms_sym <- reaction_norms %>% filter(metric == "sym.cm2")

reaction_norm_sym <- ggplot(norms_sym, aes(x=final, y=mean)) +
  # DATA 
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  # AESTHETICS 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        legend.position = "none", 
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Depth (m)", y = "Symbiont Density (cells/cm²)") +
  scale_y_continuous(labels = label_scientific(digits = 1)) 

reaction_norm_sym 

### CHL 
norms_chl <- reaction_norms %>% filter(metric == "chla.ug.cm2")

reaction_norm_chl <- ggplot(norms_chl, aes(x=final, y=mean)) +
  # DATA 
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  # AESTHETICS 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        legend.position=c(0.45, 0.85), 
        legend.background = element_rect(fill = "white", color = NA),
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_blank(), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),         
        plot.margin=unit(c(0,0,0,1),"cm")) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Depth (m)", y = "Chlorophyll (ug/cm²)", color = "Treatment") #(cm²)

reaction_norm_chl 

### FvFm
norms_fvfm <- reaction_norms %>% filter(metric == "Fv.Fm")

reaction_norm_fvfm <- ggplot(norms_fvfm, aes(x=final, y=mean)) +
  # DATA 
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  # AESTHETICS 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        legend.position = "none", 
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Depth (m)", y = "Fv/Fm") + #(cm²) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) 

reaction_norm_fvfm

### Qm 
norms_qm <- reaction_norms %>% filter(metric == "Qm")

reaction_norm_qm <- ggplot(norms_qm, aes(x=final, y=mean)) +
  # DATA 
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  # AESTHETICS 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        legend.position = "none", 
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_blank(), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),         
        plot.margin=unit(c(0,0,0,0),"cm")) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Depth (m)", y = "Qm") #(cm²)

reaction_norm_qm

### biomass
norms_biomass <- reaction_norms %>% filter(metric == "Host_AFDW_mg.cm2")

reaction_norm_biomass <- ggplot(norms_biomass, aes(x=final, y=mean)) +
  # DATA 
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  # AESTHETICS 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        legend.position = "none", 
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  labs(x = "Depth (m)", y = "Host Biomass (mg/cm²)") + #(cm²) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) 

reaction_norm_biomass

### protein
norms_prot <- reaction_norms %>% filter(metric == "prot_mg.cm2")

reaction_norm_prot <- ggplot(norms_prot, aes(x=final, y=mean)) +
  # DATA 
  geom_line(aes(group = group), size = 3) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color=home_away), width = 0.05) +
  geom_point(aes(color=home_away), size = 8) + 
  # AESTHETICS 
  facet_wrap(~group, ncol = 3, labeller = labeller(group = facet_labels))+#, strip.position = "left") + 
  scale_color_manual(values = custom_colors) +  # Set custom point colors
  theme_classic2()+
  theme(text = element_text(size=30),
        legend.position = "none", 
        strip.background = element_blank(), strip.border = element_blank(), strip.text = element_blank(), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),         
        plot.margin=unit(c(0,0,0,0),"cm")) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
  labs(x = "Depth (m)", y = "Protein (mg/cm²)")

reaction_norm_prot 

### Arrange
norms_photophys <- grid.arrange(reaction_norm_sym, reaction_norm_chl, nrow=2,
                       bottom = grobTree(
                         textGrob("*", x = 0.32, y = 136, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), # SYM 1
                         textGrob("*", x = 0.57, y = 51, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), # chl 2
                         textGrob("*", x = 0.84, y = 49, just = "left", gp = gpar(fontsize = 50, fontface = "bold")) # chl 3
                       ))

norms_photosynth <- grid.arrange(reaction_norm_fvfm, reaction_norm_qm, nrow=2,
                        bottom = grobTree(
                          textGrob("*", x = 0.27, y = 158, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          textGrob("*", x = 0.53, y = 142, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          textGrob("*", x = 0.82, y = 142, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          
                          textGrob("*", x = 0.3, y = 48, just = "left", gp = gpar(fontsize = 50, fontface = "bold")),
                          textGrob("*", x = 0.59, y = 53, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          textGrob("*", x = 0.87, y = 53, just = "left", gp = gpar(fontsize = 50, fontface = "bold")) 
                        ))

norms_tissue <- grid.arrange(reaction_norm_biomass, reaction_norm_prot, nrow=2,
                        bottom = grobTree(
                          textGrob("*", x = 0.26, y = 140, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          textGrob("*", x = 0.55, y = 138, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          textGrob("*", x = 0.85, y = 132, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          
                          textGrob("*", x = 0.26, y = 52, just = "left", gp = gpar(fontsize = 50, fontface = "bold")),
                          textGrob("*", x = 0.55, y = 60, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                          textGrob("*", x = 0.85, y = 46, just = "left", gp = gpar(fontsize = 50, fontface = "bold")) 
                        ))
           

###
norms_photophys <- grid.arrange(reaction_norm_sym, reaction_norm_chl, nrow=2,
                                bottom = grobTree(
                                  textGrob("*", x = 0.3, y = 158, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), # SYM 1
                                  textGrob("*", x = 0.56, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), # chl 2
                                  textGrob("*", x = 0.85, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")) # chl 3
                                ))

norms_photosynth <- grid.arrange(reaction_norm_fvfm, reaction_norm_qm, nrow=2,
                                 bottom = grobTree(
                                   textGrob("*", x = 0.26, y = 160, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                                   textGrob("*", x = 0.55, y = 160, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                                   textGrob("*", x = 0.85, y = 160, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                                   
                                   textGrob("*", x = 0.26, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")),
                                   textGrob("*", x = 0.55, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                                   textGrob("*", x = 0.85, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")) 
                                 ))

norms_tissue <- grid.arrange(reaction_norm_biomass, reaction_norm_prot, nrow=2,
                             bottom = grobTree(
                               textGrob("*", x = 0.26, y = 160, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                               textGrob("*", x = 0.55, y = 160, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                               textGrob("*", x = 0.85, y = 160, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                               
                               textGrob("*", x = 0.26, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")),
                               textGrob("*", x = 0.55, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")), 
                               textGrob("*", x = 0.85, y = 65, just = "left", gp = gpar(fontsize = 50, fontface = "bold")) 
                             ))

                    
norms_phys <- grid.arrange(norms_photophys, norms_photosynth, norms_tissue, nrow=1)

ggsave("TL_TRANS_reaction_norms_phys.jpg", plot = norms_phys, path = '~/Desktop/TL_TRANS_Manuscript/Graphs/', width = 25, height = 18)


# Fig Sx. Boxplots --------------------------------------------------------

# pivot longer 
full_phys_long <- full_phys %>% pivot_longer(cols = c("A", "CA", "di", "Cdi", "D2", "Dsmall", "Dadult", "chla.ug.cm2", "sym.cm2", "prot_mg.cm2", "Host_AFDW_mg.cm2", "Sym_AFDW_mg.cm2", "CM2.year", "Qm", "Fv.Fm"), names_to = "metric", values_to = "value")

# order 
full_phys_long$metric <- factor(full_phys_long$metric,levels = c("percent_survival", "CM2.year", "D2", "Dsmall", "Dadult",
             "A", "di", "CA", "Cdi",
             "Fv.Fm", "Qm", "prot_mg.cm2", "Host_AFDW_mg.cm2",
             "chla.ug.cm2", "sym.cm2", "Sym_AFDW_mg.cm2"))
full_phys_long$final <- factor(full_phys_long$final, levels = c("5", "18" )) # Replace with desired order
full_phys_long$home_away <- factor(full_phys_long$home_away, levels = c("Home", "Away" ))

# make labels 
facet_labels <- c(
  "D2" = "Corallite Density (per/cm²)", 
  "Dsmall" = "Bud Density (per/cm²)",
  "Dadult" = "Adult Density (per/cm²)",
  "CA" = "Calice Area (cm²)", 
  "Cdi" = "Calice Diameter (cm)",
  "CM2.year" = "Growth (cm²/year)",  
  "A" = "Corallite Area (cm²)", 
  "di" = "Corallite Diameter (cm)", 
  "percent_survival" = "Percent Survival (%)", 
  "chla.ug.cm2" = "Chlorophyll α (μg/cm²)", 
  "Fv.Fm" = "Fv.Fm", 
  "Qm" = "Qm", 
  "sym.cm2" = "Symbiont Density (cells/cm²)",
  "prot_mg.cm2" = "Soluble Protein (mg/cm²)", 
  "Host_AFDW_mg.cm2" = "Host Biomass (mg/cm²)",
  "Sym_AFDW_mg.cm2" = "Symbiont Biomass (mg/cm²)"   
)

### OFAV SHALLOW 

data_OFAVS <- full_phys_long %>%
  filter(group == "OFAV_S") %>%
  filter(metric != "NA")

box_OFAVS <- ggplot(data_OFAVS, aes(x=final, y=value, fill=home_away)) + 
  facet_wrap( ~ metric, scales = "free", labeller = labeller(metric = facet_labels), ncol = 8, strip.position= "left") +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  theme_classic2() + 
  theme(text = element_text(size=11),
    legend.position = "none",
    strip.placement = "outside",
    plot.title = element_text(face="bold", hjust = 0.5),
    axis.title.x = element_text(face="bold"),
    plot.margin=unit(c(0.1,0.1,0.1,0),"cm")
  ) + 
  labs(x = "Depth (m)", y = "", title = "OFAV Shallow")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",  tip.length = 0,
                     label.y.npc = 0.83, label.x.npc = 0.4, 
                     size = 11, color = "blue", 
                     symnum.args = list(cutpoints = c(0,  0.05, Inf), symbols = c("*", "")))  

### OFAV DEEP 

data_OFAVP <- full_phys_long %>%
  filter(group == "OFAV_P") %>%
  filter(metric != "NA")

box_OFAVP <- ggplot(data_OFAVP, aes(x=final, y=value, fill=home_away)) + 
  facet_wrap( ~ metric, scales = "free", labeller = labeller(metric = facet_labels), ncol = 8, strip.position= "left") +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  theme_classic2() + 
  theme(text = element_text(size=11),
        legend.position=c(0.93, 0.2), 
        legend.background = element_rect(fill = "white", color = NA),
        strip.placement = "outside",
        plot.title = element_text(face="bold", hjust = 0.5),
        axis.title.x = element_text(face="bold"),
        plot.margin=unit(c(0.5,0.1,0.1,0),"cm")
  ) + 
  labs(x = "Depth (m)", y = "", title = "OFAV Deep", fill = "Treatment")+
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",  tip.length = 0,
                     label.y.npc = 0.83, label.x.npc = 0.4, 
                     size = 11, color = "blue", 
                     symnum.args = list(cutpoints = c(0,  0.05, Inf), symbols = c("*", ""))) 
### OFRA DEEP 

data_OFRAP <- full_phys_long %>%
  filter(group == "OFRA_P") %>%
  filter(metric != "NA")

box_OFRA <- ggplot(data_OFRAP, aes(x=final, y=value, fill=home_away)) + 
  facet_wrap( ~ metric, scales = "free", labeller = labeller(metric = facet_labels), ncol = 8, strip.position= "left") +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  theme_classic2() + 
  theme(text = element_text(size=11),
        legend.position = "none",
        strip.placement = "outside",
        plot.title = element_text(face="bold", hjust = 0.5),
        axis.title.x = element_text(face="bold"),
        plot.margin=unit(c(0.5,0.1,0.1,0),"cm")
  ) + 
  labs(x = "Depth (m)", y = "", title = "OFRA Deep") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",  tip.length = 0,
                     label.y.npc = 0.83, label.x.npc = 0.4, 
                     size = 11, color = "blue", 
                     symnum.args = list(cutpoints = c(0,  0.05, Inf), symbols = c("*", "")))  

### JOIN & SAVE 

phys_arrange <- grid.arrange(
  box_OFAVS, box_OFAVP, box_OFRA,  # Use the grob versions
  nrow = 3
  #widths = c(0.63, 0.10, 0.27),
  #bottom = text_grob("Days Since Transplantation", size = 30, face = "bold")
)  

ggsave("TL_TRANS_FigSx_boxplots.jpg", plot = phys_arrange, path = '~/Desktop/TL_TRANS_Manuscript/Graphs/', width = 10, height =13)


# X Physiology graphs ---------------------------------------------------------------------

# Graph Qm
Qm <- ggplot(PAM, aes(x=factor(full_treatment, level=x_order), y=Qm, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(0.6, 0.6, 0.6)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(y= bquote("Qm"), x = "", fill='Origin') 

Qm
pairwise.wilcox.test(PAM$Qm, PAM$full_treatment, p.adjust.method="none")

# Graph Fv/Fm
fvfm <- ggplot(PAM, aes(x=factor(full_treatment, level=x_order), y=Fv.Fm, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(0.75, 0.75, 0.75)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(y= bquote("Fv/Fm"), x = bquote(""), fill='Origin') + 
  scale_x_discrete(labels=c('Shallow', 'Deep', 'Shallow', 'Deep', 'Shallow','Deep')) 

fvfm
pairwise.wilcox.test(PAM$Fv.Fm, PAM$full_treatment, p.adjust.method="none")

# Chlorophyll a 
chla <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=chla.ug.cm2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(15, 15, 15)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        #axis.text.x = element_blank(), 
        #axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "Depth Treatment", y = bquote("Chlorophyll α (μg/" ~ cm^2~ ")"), fill='Origin') + 
  scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
chla
pairwise.wilcox.test(raw$chla.ug.cm2, raw$full_treatment, p.adjust.method="none")

# Symbiont Density 
sym <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=sym.cm2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(4000000, 4000000, 4000000)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "", y = bquote("Symbiont Density (cells/" ~ cm^2~ ")"), fill='Origin')
sym
pairwise.wilcox.test(raw$sym.cm2, raw$full_treatment, p.adjust.method="none")

# Protein 
protein <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=prot_mg.cm2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(2.3, 2.3, 2.3)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(1,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "", y = bquote("Protein (mg/" ~ cm^2~ ")"), fill='Origin') 
#scale_x_discrete(labels=c('Shallow', 'Deep', 'Shallow', 'Deep', 'Shallow', 'Deep')) 
protein
pairwise.wilcox.test(raw$prot_mg.cm2, raw$full_treatment, p.adjust.method="none")

# Symbiont Biomass 
AFDW_sym <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=Sym_AFDW_mg.cm2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(40, 40, 40)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        #axis.text.x = element_blank(), 
        #axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "Depth Treatment", y = bquote("Symbiont Biomass (mg/" ~ cm^2~ ")"), fill='Origin') + 
  scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
AFDW_sym
pairwise.wilcox.test(raw$Sym_AFDW_mg.cm2, raw$full_treatment, p.adjust.method="none")

# Host Biomass 
AFDW_host <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=Host_AFDW_mg.cm2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(16, 16, 16)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "", y = bquote("Host Biomass (mg/"~cm^2~")"), fill='Origin') 

AFDW_host
pairwise.wilcox.test(raw$Host_AFDW_mg.cm2, raw$full_treatment, p.adjust.method="none")

# X Physiology arrange ------------------------------------------------------

blank <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=Host_AFDW_mg.cm2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  #geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  # AESTHETICS 
  theme_classic() +
  theme(legend.position=c(0.35, 0.3), 
        legend.background = element_rect(fill = "white", color = NA), # White background
        legend.margin = margin(200, 150, 400,150),
        text = element_text(size=40), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        legend.key.size = unit(3, "cm"),       # Increase overall size of the legend keys
        plot.margin=unit(c(1,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels = c(expression(italic("O. faveolata ") * "Shallow"), 
                                                                            expression(italic("O. faveolata ") * " Deep"),
                                                                            expression(italic("O. franksi ") * " Deep")) ) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(y= bquote(""), x = "", fill='Original Population')

phys_col1 <- plot_grid(protein, AFDW_host, 
                       Qm, fvfm,
                       sym, AFDW_sym, 
                       chla, blank,    
                       ncol = 2, align = "vh", 
                       labels = c("A","B", "C", "D", "E", "F",  "G",""),  label_size = 40, label_x = 0.04, label_y = 0.99)

#save 2x6 graphs 
ggsave("TL_TRANS_arrange_phys.jpg", plot = phys_col1, path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 18, height = 40)

# X Morphology --------------------------------------------------------------

# growth 
growth <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=CM2.year, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(6,6,6)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(1.1,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "", y = bquote("Growth ("~cm^2~"/year)"), fill='Origin')
#scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
growth
pairwise.wilcox.test(raw$CM2.year, raw$full_treatment, p.adjust.method="none")

# Corallite Density  
density <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=D2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(8,8,8)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(1.1,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "", y = bquote("Corallite Density (calices/"~cm^2~")"), fill='Origin') 
#scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
density 
pairwise.wilcox.test(raw$D2, raw$full_treatment, p.adjust.method="none")

# Corallite Density  
density_small <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=Dsmall, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(3,3,3)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        #axis.text.x = element_blank(), 
        #axis.ticks.x= element_blank(),
        plot.margin=unit(c(1.1,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "Depth Treatment", y = bquote("Bud Density (calices/"~cm^2~")"), fill='Origin') +
  scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
density_small 
pairwise.wilcox.test(raw$Dsmall, raw$full_treatment, p.adjust.method="none")

# Calice Area 
area <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=A, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(0.21, 0.21, 0.21)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "", y = bquote("Coralite Area ("~cm^2~")"), fill='Origin') 
#scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
area
pairwise.wilcox.test(raw$A, raw$full_treatment, p.adjust.method="none")

# Columella area 
Carea <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=CA, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(0.09, 0.09, 0.09)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        #axis.text.x = element_blank(), 
        #axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "Depth Treatment", y = bquote("Calice Area ("~cm^2~")"), fill='Origin') + 
  scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 

Carea
pairwise.wilcox.test(raw$CA, raw$full_treatment, p.adjust.method="none")

# Calice diameter 
diameter <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=di, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(0.58, 0.58, 0.58)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "", y = bquote("Coralite Diameter (cm)"), fill='Origin') 
#scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
diameter
pairwise.wilcox.test(raw$di, raw$full_treatment, p.adjust.method="none") 

# Columnella diameter
Cdiameter <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=Cdi, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  stat_compare_means(comparisons = treatment_comparisons, method = "wilcox.test", size = 8, label.y = c(0.34, 0.34, 0.34)) + 
  # AESTHETICS 
  theme_classic() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        #axis.text.x = element_blank(), 
        #axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels=c('OFAV Shallow', 'OFAV Deep', 'OFRA Deep')) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(x = "Depth Treatment", y = bquote("Calice Diameter (cm)"), fill='Origin') + 
  scale_x_discrete(labels=c('S', 'D', 'S', 'D', 'S', 'D')) 
Cdiameter
pairwise.wilcox.test(raw$Cdi, raw$full_treatment, p.adjust.method="none")

# X morphology arrange  -----------------------------------------------------

blank2 <- ggplot(raw, aes(x=factor(full_treatment, level=x_order), y=Host_AFDW_mg.cm2, fill=factor(group, level=measurement_order))) +
  # DATA 
  geom_boxplot(alpha=0.6, outlier.size=0) +
  #geom_point(aes(color=factor(group, level=measurement_order)), position = pj, size=2, show.legend = FALSE)+
  # AESTHETICS 
  theme_classic() +
  theme(legend.position=c(0.3, 0.6), 
        legend.background = element_rect(fill = "white", color = NA), # White background
        legend.margin = margin(30, 0, 200,80),
        text = element_text(size=40), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        legend.key.size = unit(3, "cm"),       # Increase overall size of the legend keys
        plot.margin=unit(c(1,0,0,1),"cm")
  ) + 
  scale_fill_manual(values = c("#EBB134", "#ED6B5F", "#6060A8"), labels = c(expression(italic("O. faveolata ") * "Shallow"), 
                                                                            expression(italic("O. faveolata ") * " Deep"),
                                                                            expression(italic("O. franksi ") * " Deep")) ) + 
  scale_color_manual(values = c("#EBB134", "#ED6B5F", "#6060A8")) + 
  labs(y= bquote(""), x = "", fill='Original Population')

phys_col1 <- plot_grid(protein, blank, sym, chla, Qm, fvfm, AFDW_sym, AFDW_host, 
                       ncol = 2, align = "vh", 
                       labels = c("A","", "B", "C", "D", "E", "F",  "G"),  label_size = 35, label_x = 0.18, label_y = 0.99)


morph_col1 <- plot_grid(growth, density, density_small, area, diameter, blank2, Carea, Cdiameter, blank2, 
                        ncol = 3, align = "vh", 
                        labels = c("A", "B", "C", "D", "E", "", "G", "H", ""),  label_size = 40, label_x = 0.04, label_y = 0.99)

#save 2x6 graphs 
ggsave("TL_TRANS_arrange_morph.jpg", plot = morph_col1, path = '~/Desktop/TL_Trans_Manuscript/Graphs/', width = 25, height = 30)
#

# X Phys & morph stats ------------------------------------------------------

# prep specific ones 
raw_OFAVP <- raw %>% filter(full_treatment %in% c("OFAV_PP", "OFAV_PS"))
raw_OFAVS <- raw %>% filter(full_treatment %in% c("OFAV_SP", "OFAV_SS"))
raw_OFRA <- raw %>% filter(full_treatment %in% c("OFRA_PP", "OFRA_PS"))
raw_OFAV_control <- raw %>% filter(full_treatment %in% c("OFAV_PP", "OFAV_SS"))
raw_species <- raw %>% filter(full_treatment %in% c("OFAV_PP", "OFRA_PP"))
columns_to_test <- c("A", "CA", "di", "Cdi", "D2","Dsmall","Dadult", "chla.ug.cm2", "sym.cm2", "prot_mg.cm2", "Host_AFDW_mg.cm2", "Sym_AFDW_mg.cm2", "CM2.year")  

# prep sepecific ones 
pam_OFAVP <- PAM %>% filter(full_treatment %in% c("OFAV_PP", "OFAV_PS"))
pam_OFAVS <- PAM %>% filter(full_treatment %in% c("OFAV_SP", "OFAV_SS"))
pam_OFRA <- PAM %>% filter(full_treatment %in% c("OFRA_PP", "OFRA_PS"))
pam_OFAV_control <- PAM %>% filter(full_treatment %in% c("OFAV_PP", "OFAV_SS"))
pam_species <- PAM %>% filter(full_treatment %in% c("OFAV_PP", "OFRA_PP"))
columns_to_test2 <- c("Qm", "Fv.Fm")  

apply_wilcoxon <- function(data, columns, group_col, group_label) {
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each column in the specified columns
  for (col in columns) {
    # Perform Wilcoxon test for each column
    test_result <- wilcox.test(data[[col]] ~ data[[group_col]], data = data, p.adjust.method = "none")
    
    # Create a summary for the test result
    summary <- data.frame(
      metric = col,      
      Comparison = group_label,
      w = test_result$statistic, 
      p = test_result$p.value)
    
    # Append the result to the results list
    results_list[[col]] <- summary}
  
  # Combine the results into a single data frame
  final_results <- do.call(rbind, results_list)
  
  return(final_results)
}

# apply to control vs. shade then control vs. deep 
results_OFAVP <- apply_wilcoxon(raw_OFAVP, columns_to_test, "full_treatment", "OFAVP")
results_OFAVS <- apply_wilcoxon(raw_OFAVS, columns_to_test, "full_treatment", "OFAVS")
results_OFRA <- apply_wilcoxon(raw_OFRA, columns_to_test, "full_treatment", "OFRA")
results_OFAV_control <- apply_wilcoxon(raw_OFAV_control, columns_to_test, "full_treatment", "OFAV Control")
results_species <- apply_wilcoxon(raw_species, columns_to_test, "full_treatment", "Species")

# apply to control vs. shade then control vs. deep 
results_pam_OFAVP <- apply_wilcoxon(pam_OFAVP, columns_to_test2, "full_treatment", "OFAVP")
results_pam_OFAVS <- apply_wilcoxon(pam_OFAVS, columns_to_test2, "full_treatment", "OFAVS")
results_pam_OFRA <- apply_wilcoxon(pam_OFRA, columns_to_test2, "full_treatment", "OFRA")
results_pam_OFAV_control <- apply_wilcoxon(pam_OFAV_control, columns_to_test2, "full_treatment","OFAV Control")
results_pam_species <- apply_wilcoxon(pam_species,columns_to_test2, "full_treatment", "Species")

# merge 
wilcox_2 <- merge(results_OFAVP, results_OFAVS, by = "metric", all = TRUE)
wilcox_3 <- merge(wilcox_2, results_OFRA, by ="metric", all = TRUE)
wilcox_4 <- merge(wilcox_3, results_OFAV_control, by ="metric", all = TRUE)
wilcox_5 <- merge(wilcox_4, results_species , by ="metric", all = TRUE)

# merge 
wilcox_pam_2 <- merge(results_pam_OFAVP, results_pam_OFAVS, by = "metric", all = TRUE)
wilcox_pam_3 <- merge(wilcox_pam_2, results_pam_OFRA, by ="metric", all = TRUE)
wilcox_pam_4 <- merge(wilcox_pam_3, results_pam_OFAV_control, by ="metric", all = TRUE)
wilcox_pam_5 <- merge(wilcox_pam_4, results_pam_species , by ="metric", all = TRUE)

# merge phys and pam 
combined_dataset <- bind_rows(wilcox_5, wilcox_pam_5)

# write to csv 
write.csv(combined_dataset, "~/Desktop/TL_TRANS_Manuscript/STATS/TL_TRANS_STATS_phys_wilcox.csv", row.names = FALSE)

## Means
master_long <- raw %>%
  pivot_longer(cols = c("A", "CA", "di", "Cdi", "D2","Dsmall", "chla.ug.cm2", "sym.cm2", "prot_mg.cm2", "Host_AFDW_mg.cm2", "Sym_AFDW_mg.cm2", "CM2.year"), names_to = "metric", values_to = "value") %>%
  filter(value != "NA") %>%
  select(colony_id, full_treatment, metric, value)

phys_means <- master_long %>%
  group_by(metric, full_treatment) %>%
  summarise(mean = signif(mean(value),3), SD = signif(sd(value),3), count = n()) %>%
  pivot_wider(names_from = full_treatment, values_from = c(mean, SD, count)) %>%
  select(metric, mean_OFAV_PP, SD_OFAV_PP, count_OFAV_PP, 
         mean_OFAV_PS, SD_OFAV_PS, count_OFAV_PS, 
         mean_OFAV_SP, SD_OFAV_SP, count_OFAV_SP,
         mean_OFAV_SS, SD_OFAV_SS, count_OFAV_SS,
         mean_OFRA_PP, SD_OFRA_PP, count_OFRA_PP,
         mean_OFRA_PS, SD_OFRA_PS, count_OFRA_PS)

master_long_pam <- PAM %>%
  pivot_longer(cols = c("Qm", "Fv.Fm"), names_to = "metric", values_to = "value") %>%
  filter(value != "NA") %>%
  select(full_treatment, metric, value)

pam_means <- master_long_pam %>%
  group_by(metric, full_treatment) %>%
  summarise(mean = signif(mean(value),3), SD = signif(sd(value),3)) %>%
  pivot_wider(names_from = full_treatment, values_from = c(mean, SD)) %>%
  select(metric, mean_OFAV_PP, SD_OFAV_PP,  mean_OFAV_PS, SD_OFAV_PS, mean_OFAV_SP, SD_OFAV_SP, mean_OFAV_SS, SD_OFAV_SS, mean_OFRA_PP, SD_OFRA_PP, mean_OFRA_PS, SD_OFRA_PS)

# merge 
physiology_means <- bind_rows(phys_means, pam_means)
write.csv(physiology_means, "~/Desktop/TL_TRANS_Manuscript/STATS/TL_TRANS_STATS_phys_means.csv", row.names = FALSE)

## Medians
phys_medians <- master_long %>%
  group_by(metric, full_treatment) %>%
  summarise(median = signif(median(value),3), count = n()) %>%
  pivot_wider(names_from = full_treatment, values_from = c(median, count)) %>%
  select(metric, 
         median_OFAV_PP, count_OFAV_PP, 
         median_OFAV_PS,  count_OFAV_PS, 
         median_OFAV_SP, count_OFAV_SP,
         median_OFAV_SS, count_OFAV_SS,
         median_OFRA_PP,  count_OFRA_PP,
         median_OFRA_PS,count_OFRA_PS)

pam_medians <- master_long_pam %>%
  group_by(metric, full_treatment) %>%
  summarise(median = signif(mean(value),3), count = n()) %>%
  pivot_wider(names_from = full_treatment, values_from = c(median, count)) %>%
  select(metric, 
         median_OFAV_PP, count_OFAV_PP, 
         median_OFAV_PS,  count_OFAV_PS, 
         median_OFAV_SP, count_OFAV_SP,
         median_OFAV_SS, count_OFAV_SS,
         median_OFRA_PP,  count_OFRA_PP,
         median_OFRA_PS,count_OFRA_PS)

# merge 
physiology_medians <- bind_rows(phys_medians, pam_medians)
write.csv(physiology_medians, "~/Desktop/TL_TRANS_Manuscript/STATS/TL_TRANS_STATS_phys_medians.csv", row.names = FALSE)

# calculate PImds 
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2745.2006.01176.x
# (Maximum median- minimum median)/ maximum median

Pmedians <- physiology_medians %>% select(metric, median_OFAV_PP, median_OFAV_PS, median_OFAV_SP, median_OFAV_SS, median_OFRA_PP, median_OFRA_PS)
PImd <- Pmedians %>%
  mutate(OFAV_P_max = max(median_OFAV_PP, median_OFAV_PS)) %>%
  mutate(OFAV_P_min = min(median_OFAV_PP, median_OFAV_PS)) %>%
  mutate(OFAV_S_max = max(median_OFAV_SP, median_OFAV_SS)) %>%
  mutate(OFAV_S_min = min(median_OFAV_SP, median_OFAV_SS)) %>%
  mutate(OFRA_P_max = max(median_OFRA_PP, median_OFRA_PS)) %>%
  mutate(OFRA_P_min = min(median_OFRA_PP, median_OFRA_PS)) %>%
  mutate(PI_OFAV_P = (OFAV_P_max - OFAV_P_min)/OFAV_P_max) %>%
  mutate(PI_OFAV_S = (OFAV_S_max - OFAV_S_min)/OFAV_S_max) %>%
  mutate(PI_OFRA_P = (OFRA_P_max - OFRA_P_min)/OFRA_P_max)

PImd2 <- PImd %>%
  select(metric, PI_OFAV_P, PI_OFAV_S,PI_OFRA_P)


PImd2$Rank <- apply(PImd2[, -1], 1, function(row) {
  values <- sort(row, decreasing = TRUE)  # Sort the row in descending order
  labels <- c("OFAV_P", "OFAV_S", "OFRA_P") 
  
  # Create a string representing the order
  order_str <- paste(labels[order(row, decreasing = TRUE)], collapse = ">")
  return(order_str)
})

idk <- PImd2 %>%
  pivot_longer(cols = c(PI_OFAV_P, PI_OFAV_S, PI_OFRA_P), names_to = "group", values_to = "PI") %>%
  group_by(group) %>%
  summarise(mean = mean(PI), 
            sum = sum(PI))

pivot_longer(cols = c("A", "CA", "di", "Cdi", "D2", "Dsmall", "chla.ug.cm2", "sym.cm2", "prot_mg.cm2", "Host_AFDW_mg.cm2", "Sym_AFDW_mg.cm2", "CM2.year", "Qm", "Fv.Fm"), names_to = "metric", values_to = "value") %>%
  
  
# X PCAs  -------------------------------------------------------------------

raw_complete <- raw %>% na.omit()
pca1 <- prcomp(raw_complete[,c(7:13, 15:20)],center=TRUE,scale.=TRUE)
pca1 <- prcomp(raw_complete[,c(13, 16:19)],center=TRUE,scale.=TRUE) # phys only no qm fvfm
pca1 <- prcomp(raw_complete[,c(7:12, 20)],center=TRUE,scale.=TRUE) # morph only 
autoplot(pca1, data=raw_complete, color = 'group', frame.type = 't', loadings=TRUE, loadings.label=TRUE) 

autoplot(
  pca1, 
  data = raw_complete, 
  colour = 'group',         # Color aesthetic for points
  frame.type = 't',         # Type of frame (e.g., ellipse)
  loadings = TRUE,          # Show PCA loadings
  loadings.label = TRUE     # Label the loadings
) +
  scale_color_manual(       # Custom colors for points
    values = c("OFAV_P" = "blue", "OFAV_S" = "red", "OFRA_P" = "green")
  ) +
  scale_fill_manual(        # Custom colors for ellipses
    values = c("OFAV_P" = "blue", "OFAV_S" = "red", "OFRA_P" = "green")
  )

# X Reaction norms ----------------------------------------------------------

# combine RAW and PAM values 
PAM <- PAM %>% mutate(colony_id = as.character(colony_id))
full_phys <- full_join(raw, PAM)
full_phys <- full_phys %>%
  select(species, group, full_treatment, final, home_away, CM2.year, D2, Dsmall,Dadult, A, CA, di, Cdi, Qm, Fv.Fm, prot_mg.cm2, sym.cm2, Sym_AFDW_mg.cm2, Host_AFDW_mg.cm2, chla.ug.cm2)

phys_longer <- full_phys %>%
  pivot_longer(cols = c("A", "CA", "di", "Cdi", "D2", "Dsmall", "Dadult", "chla.ug.cm2", "sym.cm2", "prot_mg.cm2", "Host_AFDW_mg.cm2", "Sym_AFDW_mg.cm2", "CM2.year", "Qm", "Fv.Fm"), names_to = "metric", values_to = "value") %>%
  filter(value != "NA") %>%
  select(group, final, home_away, metric, value)

reaction_norms <- phys_longer %>%
  group_by(group, final, metric, home_away)%>%
  summarise(mean = mean(value),
            sd = sd(value))

# add survival in here 
n_survival2 <- n_survival %>% 
  mutate(metric = "percent_survival") %>%
  select(final, group, metric, mean) 

reaction_norms <-  bind_rows(reaction_norms, n_survival2)

reaction_norms$metric <- factor(
  reaction_norms$metric,
  levels = c("percent_survival", "CM2.year", "D2", "Dsmall",
             "A", "di", "CA", "Cdi",
             "Fv.Fm", "Qm", "prot_mg.cm2", "Host_AFDW_mg.cm2",
             "chla.ug.cm2", "sym.cm2", "Sym_AFDW_mg.cm2"
  )) # Replace with desired order

reaction_norms$final <- factor(
  reaction_norms$final,
  levels = c("5", "18" )) # Replace with desired order

facet_labels <- c(
  "D2" = "Corallite Density (calices/cm²)", 
  "Dsmall" = "Bud Density (calices/cm²)",
  "CA" = "Calice Area (cm²)", 
  "Cdi" = "Calice Diameter (cm)",
  "CM2.year" = "Growth (cm²/year)",  
  "A" = "Corallite Area (cm²)", 
  "di" = "Corallite Diameter (cm)", 
  "percent_survival" = "Percent Survival (%)", 
  "chla.ug.cm2" = "Chlorophyll α (μg/cm²)", 
  "Fv.Fm" = "Fv.Fm", 
  "Qm" = "Qm", 
  "sym.cm2" = "Symbiont Density  (cells/cm²)",
  "prot_mg.cm2" = "Soluble Protein (mg/cm²)", 
  "Host_AFDW_mg.cm2" = "Host Biomass (mg/cm²)",
  "Sym_AFDW_mg.cm2" = "Symbiont Biomass (mg/cm²)"   
)

reaction_norm_plot <- ggplot(reaction_norms, aes(x=final, y=mean, color=group)) +
  geom_line(aes(group = group), size = 3) + 
  geom_point(size = 6) + 
  facet_wrap(~metric, scales="free", ncol = 4, labeller = labeller(metric = facet_labels), strip.position = "left") + 
  theme_bw()+
  theme(strip.placement = "outside", 
        text = element_text(size=30),
        strip.background = element_blank(), 
        strip.border = element_blank(), 
        strip.text = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = c(0.88, 0.08) ) + 
  scale_color_manual(values = c("#ED6B5F","#EBB134", "#6060A8") ,labels = c(expression(italic("O. faveolata ") * " Deep"), 
                                                                            expression(italic("O. faveolata ") * "Shallow"), 
                                                                            expression(italic("O. franksi ") * " Deep")) ) + 
  labs(x = "Treatment", y = "", color = "Original Population")
reaction_norm_plot

ggsave("TL_TRANS_reaction_norms.jpg", plot = reaction_norm_plot, path = '~/Desktop/TL_TRANS_Manuscript/Graphs/', width = 20, height = 26)

### MORPH NORMS

# get values for means 
x_cal_a <- mean(raw$CA, na.rm = TRUE)
x_cal_di <- mean(raw$Cdi, na.rm = TRUE)
x_cor_a <- mean(raw$A, na.rm = TRUE)
x_cor_di <- mean(raw$di, na.rm = TRUE)

# standardize 
morph_reaction_norms <- reaction_norms %>%
  filter(metric %in% c("di", "Cdi", "A", "CA")) %>%
  mutate(std_value = case_when(
    metric == "di"  ~ mean - x_cor_di,
    metric == "Cdi" ~ mean - x_cal_di,
    metric == "A"   ~ mean - x_cor_a,
    TRUE            ~ mean - x_cal_a  # This acts as the "else" case
  ))

ggplot(morph_reaction_norms, aes(x=final, y = std_value, color = home_away)) +
  geom_point() + 
  facet_wrap( ~ metric + group, nrow=4)
#add error bars 


min_val <- min(raw$A, na.rm = TRUE)
max_val <- max(raw$A, na.rm = TRUE)
raw$CA_norm <- (raw$A - min_val)/(max_val - min_val)

norms <- raw %>%
  group_by(group, home_away) %>%
  summarise(mean = mean(CA_norm, na.rm = TRUE),
            sd = sd(CA_norm, na.rm = TRUE))
norms <- norms %>%
  pivot_wider(names_from = home_away, values_from = c(mean, sd)) %>%
  mutate(mean_Away = mean_Away - mean_Home) 

#%>%
mutate(mean_Home = 0) 

norms <- norms %>% 
  pivot_longer(-group, 
               names_to = c(".value", "home_away"), 
               names_sep="_" )


norms$home_away <- factor(norms$home_away, levels = c("Home", "Away"))

ggplot(norms, aes(x=home_away, y=mean, color=group, group=group))+ 
  geom_point(size=3) + 
  geom_line(size=2) +
  facet_wrap(~group) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.05) + 
  theme_classic2()

