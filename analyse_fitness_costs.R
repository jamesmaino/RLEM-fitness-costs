library(tidyverse)
library(lme4)
library(broom)
library(broom.mixed)
library(lmerTest)
library(AICcmodavg)

# to do 
# reanalysis
# compare ddpcr with amplicon

#issues
# did not match the existing data

# Load and clean data.
# online spread sheet
# d0 = read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQmErDprvdoLAIFM9Lf9Fr7JnH3_2nOYlBmYwTwr814uPDgecxpgWnExdj4WCuCVVDd915GskIWf_NY/pub?gid=0&single=true&output=csv")

# archived spreadsheet
d0 = "data/RLEM database for M7, M8 fitness and RMS UM00057.xlsx" %>%
  readxl::read_xlsx(sheet="data")

amp = readxl::read_xlsx("data/Ant_RLEM_MiSeq_v3_svs_JM.xlsx", sheet="TRANSPOSED") %>% 
  rename(location = "Field Site") %>% 
  mutate(date = as.Date(Date, format="%d/%m/%Y")) %>%
  select(location, date,  plot=Plot, quadrat=Replicate, R) %>% 
  mutate(quadrat = as.integer(quadrat))

controls = d0 %>%
  mutate(trtcode = ifelse(is.na(trtcode), 1, trtcode)) %>%
  filter(trtcode == 1) %>% 
  distinct(plot) %>%
  pull(plot)

d1  = d0 %>% 
  # samples were taken over multiple dates so rough date is used which is the last date of sampling for a given field trip
  mutate(rough_date = as.Date(rough_date)) %>%
  
  # get mean across quadrats where multiple samples are available
  dplyr::select(date=rough_date,	location,	month, plot,	
                trtcode, damage, contains("kdr")) %>%
  gather(var, val, -date,	-location,	-month,	-plot,	-trtcode,	-damage) %>% 
  filter(var != "kdr_percent.mean") %>%
  separate(var, c('var', 'quadrat'), sep = '\\.q', convert = TRUE) %>% 
  filter(!is.na(quadrat)) %>%
  filter(!is.na(val)) %>%
  filter(var %in% c('kdr_percent', 'kdr_insects.tested')) %>% 
  spread(var, val, convert = TRUE) %>% 
  mutate(kdr_percent = as.numeric(kdr_percent)) %>%
  mutate(kdr_prop_ddpcr = kdr_percent/100) %>%
  
  # join to amplicon reanalysis
  full_join(amp) %>%
  mutate(kdr_prop = R) %>%
  # assume 50 if missing
  # mutate(kdr_insects.tested = 
  #          ifelse(is.na(kdr_insects.tested), 50, kdr_insects.tested)) %>%
  mutate(alleles_screened = 2*kdr_insects.tested) %>%
  # filter(!is.na(kdr_percent)) %>% 
  mutate(year = format(date, '%Y')) %>% 
  mutate(years = as.integer(year)- 2017) %>% 
  mutate(days = as.numeric(date - as.Date(ifelse(location == 'Arthur River', '2017-08-01', '2017-06-26')))) %>% 
  mutate(days_s = (days - mean(days))/sd(days)) %>% 
  
  # remove treatment plots from tintinara
  filter(plot %in% controls)


# there are mismatches in the data. 
d1 %>% 
  select(location, date,  plot, quadrat, 
         kdr_insects.tested , kdr_prop, kdr_prop_ddpcr) %>%
  filter(is.na(kdr_prop)) %>% 
  # filter(location=="Arthur River") %>%
  View

# plot correlation between two kdr estimation methods
d1 %>%
  filter(location == "Tintinara") %>%
  ggplot(aes(kdr_prop_ddpcr, kdr_prop)) +
  geom_point(alpha=0.5) + 
  ylab("Amplicon R frequency estimate") + 
  xlab("ddPCR R frequency estimate") + 
  geom_abline(intercept = 0, slope=1) + 
  theme_bw()
ggsave("plots/explore_amplicon/ddPCRcorrelation_Tintinara.png", height=4, width=5)

# r2
lm(kdr_prop_ddpcr ~ kdr_prop, filter(d1, location == "Tintinara") ) %>% summary

# show top outliers
d1 %>% 
  drop_na(kdr_prop_ddpcr, kdr_prop) %>%
  mutate(resid = lm(kdr_prop_ddpcr ~ kdr_prop, d1)$residuals) %>% 
  arrange(desc(resid)) %>%
  write_csv("plots/explore_amplicon/descending_outliers.csv")




# Fit a logistic regression with random effects for repeated measured on plot. kdr proportion is the response variable weighted by the number of alleles screened (individual mites X2).  

d = d1 %>% 
  drop_na(kdr_prop, alleles_screened)

# Full model includes fixed effect of time and location (with interaction) and random effect of plot. 
d %>% distinct(location, plot) %>% View

m_null = 
  glmer(kdr_prop~  (1|plot), 
        weights = alleles_screened, data=d, family = binomial)
m_loc = 
  glmer(kdr_prop~ location  + (1|plot), 
        weights = alleles_screened, data=d, family = binomial)
m_time = 
  glmer(kdr_prop~ days_s  + (1|plot), 
        weights = alleles_screened, data=d, family = binomial)
m_loc_time = 
  glmer(kdr_prop~ location + days_s  + (1|plot),
        weights = alleles_screened, data=d, family = binomial)
m_loc_time_int = 
  glmer(kdr_prop~ location/days_s -1  + (1|plot), 
        weights = alleles_screened, data=d, family = binomial)

# model comparison
AICtable = aictab(list(
  `Null model` = m_null, 
  `Loc`        = m_loc,
  `Time`       = m_time,
  `Loc + Time `= m_loc_time,
  `Loc + Time + Loc*Time `= m_loc_time_int
))
AICtable

write_csv(AICtable, "plots/model_comparison.csv")

# output model coefficients
summary(m_loc_time_int)
invlogit = function(x) exp(x)/(1+exp(x))
m1coeff = broom::tidy(m_loc_time_int, conf.int=TRUE) %>% 
  mutate(MS_string = 
           sprintf("%1.2f [%1.2f. %1.2f]", estimate, conf.low, conf.high)) %>%
  # mutate(estimate.r = invlogit(estimate)*100, 
  #        conf.low.r = invlogit(conf.low)*100, 
  #        conf.high.r = invlogit(conf.high)*100) %>%
  # mutate(MS_string.real = 
  #           sprintf("%1.2f [%1.2f. %1.2f]", estimate.r, conf.low.r, conf.high.r)) %>%
  identity
write_csv(m1coeff, "plots/model_coefficients.csv")


# Ftest time
anova(m_time, m_null, test="Chisq")

# Ftest location
anova(m_loc, m_null, test="Chisq")

# predict data
d$pred_m1 = predict(m_loc_time_int, type= "response") # predictions with random effects
d$pred_fixed_m1 = predict(m_loc_time_int, type= "response", re.form=~0) # fixed effect only


# Comparison of AIC scores and liklihood ratio tests, suggests an effect of time and location of resistance frequencies but no interation of location and decrease in resistance through time i.e. no evidence for different fitness costs by location.   


# some diagnostic plots
plot(m_loc_time_int, type = c("p", "smooth"))
plot(m_loc_time_int, sqrt(abs(resid(.))) ~ fitted(.),
     type = c("p", "smooth"))
lattice::qqmath(m_loc_time_int, id = 0.05) 

# Summarise results across all replicates and plot means and standard error. See supplementary plot for individual replicates and random effects. 
dsum = d %>%
  group_by(date, location, trtcode) %>% 
  summarise(
    days = mean(days),
    days_s = mean(days_s),
    kdr_se = sd(kdr_percent, na.rm=T)/sqrt(n()), 
    kdr_percent = mean(kdr_percent, na.rm=T), 
    n = n(), .groups="drop"
  )
dsum$pred_fixed_m1 =  predict(m_loc_time_int, type="response", re.form=~0, newdata = dsum)

p2 = 
  dsum %>%
  mutate(`Trial site` = location) %>%
  ggplot(aes(x=date, y=kdr_percent, colour=`Trial site`)) +
  geom_point(aes(shape=`Trial site`)) +
  geom_line(aes(y=100*pred_fixed_m1, group = `Trial site`), size=1) + 
  geom_errorbar(aes(ymin=kdr_percent-kdr_se, ymax=kdr_percent+kdr_se),
                width = 10) + 
  xlab('Date') + 
  ylab(expression("Average "*italic(kdr)*' allele frequency (%)')) +
  theme_bw() + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "6 months") +
  scale_color_viridis_d(begin=0.3, end = 0.8, alpha=0.7)
p2
ggsave(plot=p2, 'plots/Fig3_kdr_summary1.jpeg', width=6, height=4)

# From here we can fit a population genetics function for expected allele frequency through time. Using this function we can estimate the parameter for relative fitness of resistance genotypes. Starting resistance frequencies for the two locations are left as model parameters due to uncertainty in measurements.   
# 
# First we need to convert sample time to generation number, assuming generation 1 was the first sample point.

# assume first measure was generation 1
gens_per_year = 3  # assume 3 generations per year and linear development
diapause_days = 181  # assume no development from Nov-April

d = d %>%
  mutate(
    generation = 1 + 
      round(
        gens_per_year * (days - years*diapause_days)/
          (365 - diapause_days)))  %>% 
  group_by(location, plot) %>% 
  mutate(initR = mean(kdr_prop[days <= 3], na.rm=T))

# Next we need a function that can estimate resistance frequency given the starting resistance level, number of generations, dominance of fitness costs, and fitness costs. The function is based on some [simple pop gen equations]( https://www.radford.edu/~rsheehy/Gen_flash/ABLE_Workshop/Popgen_Equations.pdf):
  
estimate_R_freq = function(w_RR, h = 0.5, generation, location, initR_A, initR_T){
  # w_RR: relative fitness of RR genotype
  #h: dominance
  w_SS = 1 # relative fitness of SS genotype
  w_RS = w_RR * h + (1 - h)*w_SS # relative fitness of RS genotype 
  
  n = length(generation)
  p_out = numeric(n)
  for(i in 1:n) {
    p =  ifelse(location[i]=="Tintinara", initR_T, initR_A)# intial R frequency
    current_gen = 1
    while(current_gen < generation[i]){
      q = 1 - p # S frequency
      w_bar = p^2*(w_RR) + 2*p*q*(w_RS) + q^2*(w_SS) # mean fitness
      p = p^2 * w_RR / w_bar + p * q * w_RS / w_bar # update R frequency
      current_gen = current_gen + 1
    }
    p_out[i] = p
  }
  return(p_out)
}
estimate_R_freq(w_RR=0.7, h=0, generation=2,  
                location="Tintinara", initR_A=0.9, initR_T=0.5)

# We can then leave fitness cost as a free parameter and estimate it from the data. I fit three models assuming fitness is recessive  (h = 0), intermediate recessive (h = 0.5), and dominant (h = 1).
# 

m2_recessive = 
  nls(kdr_prop~estimate_R_freq(w_RR, h=0, generation,  location, initR_A, initR_T), 
      data=d, start=c(w_RR=0.5, initR_A=0.9, initR_T=0.5))
m2_intermediate = 
  nls(kdr_prop~estimate_R_freq(w_RR, h=0.5, generation,  location, initR_A, initR_T), 
      data=d, start=c(w_RR=0.5, initR_A=0.9, initR_T=0.5))
m2_dominant = 
  nls(kdr_prop~estimate_R_freq(w_RR, h=1, generation,  location, initR_A, initR_T), 
      data=d, start=c(w_RR=0.5, initR_A=0.9, initR_T=0.5))

m2_free_h = 
  nls(kdr_prop~estimate_R_freq(w_RR, h, generation,  location, initR_A, initR_T), 
      data=d, start=c(w_RR=0.5, h =0.5, initR_A=0.9, initR_T=0.5))

calc_r2 = function(m){
  SS = sum((d$kdr_prop - mean(d$kdr_prop))^2)
  SSresid = sum(residuals(m)^2)
  r2 = 1 - SSresid/SS
  print(sprintf("r-squarred: %1.3f", r2))
}

summary(m2_recessive); calc_r2(m2_recessive)
summary(m2_free_h); calc_r2(m2_free_h)
summary(m2_dominant); calc_r2(m2_dominant)


broom::tidy(m2_free_h, conf.int=TRUE) %>% 
  mutate(MS_string = 
           sprintf("%1.2f [%1.2f. %1.2f]", estimate, conf.low, conf.high)) %>%
  write_csv("plots/popgen_model_coefficients.csv")


# Intermediate recessiveness marginally explains the most variance, but the models all perform similarly well. Interestingly, the estimated relative fitness of the RR genotype is fairly insensitive to dominance as shown in the above model summary outputs where relative fitness of RR varies from 0.86 to 0.89. 

# We can plot the the theoretical predictions against the data (under difference fitness dominace sceience). 

d$`h = 0.0` = predict(m2_recessive)
d$`h = 0.4` = predict(m2_free_h)
d$`h = 1.0` = predict(m2_dominant)

ggplot(d) + 
  geom_point(aes(date, kdr_percent, color = location))

dsum = d %>%
  group_by(generation, location) %>% 
  summarise(
    days = mean(days),
    days_s = mean(days_s),
    kdr_se = sd(kdr_prop, na.rm=T)/sqrt(n()), 
    kdr_prop = mean(kdr_prop, na.rm=T),
    n = n(), .groups="drop"
  ) 

pred = d %>% 
  ungroup %>%
  dplyr::select(generation, location, `h = 0.0`, 
                `h = 0.4`, `h = 1.0`) %>% distinct() %>% 
  gather(dominance, kdr_prop, `h = 0.0`, `h = 0.4`, `h = 1.0`) %>% 
  mutate(Dominance = dominance) %>% 
  mutate(`Trial site` = location)
  
p2 = dsum %>%
  mutate(`Trial site` = location) %>%
  ggplot(aes(x=generation, y=kdr_prop*100, colour=`Trial site`)) +
  geom_point(aes(shape=`Trial site`)) +
  geom_line(data = pred, aes(y=kdr_prop*100, linetype = Dominance),
            size=1) + 
  geom_errorbar(aes(ymin=(kdr_prop-kdr_se)*100, ymax=(kdr_prop+kdr_se)*100),
                width = 0.1) + 
  xlab('Generation') + 
  ylab(expression("Average "*italic(kdr)*' allele frequency (%)')) +
  theme_bw() + 
  scale_linetype_manual(values=c(1,3,5)) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_viridis_d(begin=0.3, end = 0.8, alpha=0.7) + 
  theme(legend.key.width = unit(1, "cm"))
p2
ggsave(plot=p2, 'plots/Fig4_kdr_summary2.jpg', width=6, height=4)



# plot residuals
plot(m2_recessive)
plot(m2_intermediate)
plot(m2_dominant)




# Supplementary

# We can plot kdr frequency for each replicate plot. Lines show predicted values for each plot while shapes show observed data. Size of shapes denote the number of alleles screened per sample. 

p1 = ggplot(d, 
            aes(x=days, colour=plot)) +
  # geom_text(aes(y=kdr_percent.mean, label = plot)) +
  geom_line(aes(y=pred_m1)) +
  # geom_line(aes(y=pred_fixed_m1, linetype = location), size =2, colour ='black') +
  geom_point(aes(y=kdr_percent/100, group = plot, size = alleles_screened, shape = location), 
             position = position_jitter(width = 10)) +
  xlab('days from first measurement') + 
  ylab('kdr proportion') +
  theme_classic()
p1
ggsave(plot=p1, 'plots/full_model.png', height=7, width=9)




